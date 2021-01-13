// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#include <pbrt/cpu/aggregates.h>

#include <pbrt/interaction.h>
#include <pbrt/paramdict.h>
#include <pbrt/shapes.h>
#include <pbrt/util/bits.h>
#include <pbrt/util/error.h>
#include <pbrt/util/log.h>
#include <pbrt/util/memory.h>
#include <pbrt/util/parallel.h>
#include <pbrt/util/print.h>
#include <pbrt/util/stats.h>

#include <algorithm>
#include <tuple>

namespace pbrt {

STAT_MEMORY_COUNTER("Memory/BVH", treeBytes);
STAT_RATIO("BVH/Primitives per leaf node", totalPrimitives, totalLeafNodes);
STAT_COUNTER("BVH/Interior nodes", interiorNodes);
STAT_COUNTER("BVH/Leaf nodes", leafNodes);
STAT_PIXEL_COUNTER("BVH/Nodes visited", bvhNodesVisited);

// MortonPrimitive Definition
struct MortonPrimitive {
    int primitiveIndex;
    uint32_t mortonCode;
};

// LBVHTreelet Definition
struct LBVHTreelet {
    int startIndex, nPrimitives;
    BVHBuildNode *buildNodes;
};

// BVHAggregate Utility Functions
static void RadixSort(std::vector<MortonPrimitive> *v) {
    std::vector<MortonPrimitive> tempVector(v->size());
    constexpr int bitsPerPass = 6;
    constexpr int nBits = 30;
    static_assert((nBits % bitsPerPass) == 0,
                  "Radix sort bitsPerPass must evenly divide nBits");
    constexpr int nPasses = nBits / bitsPerPass;
    for (int pass = 0; pass < nPasses; ++pass) {
        // Perform one pass of radix sort, sorting _bitsPerPass_ bits
        int lowBit = pass * bitsPerPass;
        // Set in and out vector pointers for radix sort pass
        std::vector<MortonPrimitive> &in = (pass & 1) ? tempVector : *v;
        std::vector<MortonPrimitive> &out = (pass & 1) ? *v : tempVector;

        // Count number of zero bits in array for current radix sort bit
        constexpr int nBuckets = 1 << bitsPerPass;
        int bucketCount[nBuckets] = {0};
        constexpr int bitMask = (1 << bitsPerPass) - 1;
        for (const MortonPrimitive &mp : in) {
            int bucket = (mp.mortonCode >> lowBit) & bitMask;
            CHECK_GE(bucket, 0);
            CHECK_LT(bucket, nBuckets);
            ++bucketCount[bucket];
        }

        // Compute starting index in output array for each bucket
        int outIndex[nBuckets];
        outIndex[0] = 0;
        for (int i = 1; i < nBuckets; ++i)
            outIndex[i] = outIndex[i - 1] + bucketCount[i - 1];

        // Store sorted values in output array
        for (const MortonPrimitive &mp : in) {
            int bucket = (mp.mortonCode >> lowBit) & bitMask;
            out[outIndex[bucket]++] = mp;
        }
    }
    // Copy final result from _tempVector_, if needed
    if (nPasses & 1)
        std::swap(*v, tempVector);
}

// BVHSplitBucket Definition
struct BVHSplitBucket {
    int count = 0;
    Bounds3f bounds;
};

// BVHPrimitive Definition
struct BVHPrimitive {
    BVHPrimitive() {}
    BVHPrimitive(size_t primitiveIndex, const Bounds3f &bounds)
        : primitiveIndex(primitiveIndex),
          bounds(bounds),
          centroid(.5f * bounds.pMin + .5f * bounds.pMax) {}
    size_t primitiveIndex;
    Bounds3f bounds;
    Point3f centroid;
};

// BVHBuildNode Definition
struct BVHBuildNode {
    // BVHBuildNode Public Methods
    void InitLeaf(int first, int n, const Bounds3f &b) {
        firstPrimOffset = first;
        nPrimitives = n;
        bounds = b;
        children[0] = children[1] = nullptr;
        ++leafNodes;
        ++totalLeafNodes;
        totalPrimitives += n;
    }

    void InitInterior(int axis, BVHBuildNode *c0, BVHBuildNode *c1) {
        children[0] = c0;
        children[1] = c1;
        bounds = Union(c0->bounds, c1->bounds);
        splitAxis = axis;
        nPrimitives = 0;
        ++interiorNodes;
    }

    Bounds3f bounds;
    BVHBuildNode *children[2];
    int splitAxis, firstPrimOffset, nPrimitives;
};

// LinearBVHNode Definition
struct alignas(32) LinearBVHNode {
    Bounds3f bounds;
    union {
        int primitivesOffset;   // leaf
        int secondChildOffset;  // interior
    };
    uint16_t nPrimitives;  // 0 -> interior node
    uint8_t axis;          // interior node: xyz
};

// BVHAggregate Method Definitions
BVHAggregate::BVHAggregate(std::vector<PrimitiveHandle> p, int maxPrimsInNode,
                           SplitMethod splitMethod)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)),
      primitives(std::move(p)),
      splitMethod(splitMethod) {
    CHECK(!primitives.empty());
    // Build BVH from _primitives_
    // Initialize _bvhPrimitives_ array for primitives
    std::vector<BVHPrimitive> bvhPrimitives(primitives.size());
    for (size_t i = 0; i < primitives.size(); ++i)
        bvhPrimitives[i] = BVHPrimitive(i, primitives[i].Bounds());

    // Build BVH for primitives using _bvhPrimitives_
    // Declare _Allocator_s used for BVH construction
    pstd::pmr::monotonic_buffer_resource resource;
    Allocator alloc(&resource);
    int nThreads = MaxThreadIndex();
    std::vector<pstd::pmr::monotonic_buffer_resource> threadResources(nThreads);
    std::vector<Allocator> threadAllocators;
    for (size_t i = 0; i < nThreads; ++i)
        threadAllocators.push_back(Allocator(&threadResources[i]));

    std::atomic<int> totalNodes{0};
    std::vector<PrimitiveHandle> orderedPrims(primitives.size());
    BVHBuildNode *root;
    if (splitMethod == SplitMethod::HLBVH) {
        root = buildHLBVH(alloc, bvhPrimitives, &totalNodes, orderedPrims);
    } else {
        std::atomic<int> orderedPrimsOffset{0};
        root = buildRecursive(threadAllocators, bvhPrimitives, 0, primitives.size(),
                              &totalNodes, &orderedPrimsOffset, orderedPrims);
        CHECK_EQ(orderedPrimsOffset.load(), orderedPrims.size());
    }
    primitives.swap(orderedPrims);
    bvhPrimitives.resize(0);
    LOG_VERBOSE("BVH created with %d nodes for %d primitives (%.2f MB)",
                totalNodes.load(), (int)primitives.size(),
                float(totalNodes.load() * sizeof(LinearBVHNode)) / (1024.f * 1024.f));

    // Flatten BVH into _nodes_ array
    treeBytes += totalNodes * sizeof(LinearBVHNode) + sizeof(*this) +
                 primitives.size() * sizeof(primitives[0]);
    nodes = new LinearBVHNode[totalNodes];
    int offset = 0;
    flattenBVHTree(root, &offset);
    CHECK_EQ(totalNodes.load(), offset);
}

BVHBuildNode *BVHAggregate::buildRecursive(std::vector<Allocator> &threadAllocators,
                                           std::vector<BVHPrimitive> &bvhPrimitives,
                                           int start, int end,
                                           std::atomic<int> *totalNodes,
                                           std::atomic<int> *orderedPrimsOffset,
                                           std::vector<PrimitiveHandle> &orderedPrims) {
    DCHECK_NE(start, end);
    Allocator alloc = threadAllocators[ThreadIndex];
    BVHBuildNode *node = alloc.new_object<BVHBuildNode>();
    // Initialize _BVHBuildNode_ for primitive range
    ++*totalNodes;
    // Compute bounds of all primitives in BVH node
    Bounds3f bounds;
    for (int i = start; i < end; ++i)
        bounds = Union(bounds, bvhPrimitives[i].bounds);

    int nPrimitives = end - start;
    if (bounds.SurfaceArea() == 0 || nPrimitives == 1) {
        // Create leaf _BVHBuildNode_
        int firstPrimOffset = orderedPrimsOffset->fetch_add(nPrimitives);
        for (int i = start; i < end; ++i) {
            int index = bvhPrimitives[i].primitiveIndex;
            orderedPrims[firstPrimOffset + i - start] = primitives[index];
        }
        node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
        return node;

    } else {
        // Compute bound of primitive centroids, choose split dimension _dim_
        Bounds3f centroidBounds;
        for (int i = start; i < end; ++i)
            centroidBounds = Union(centroidBounds, bvhPrimitives[i].centroid);
        int dim = centroidBounds.MaxDimension();

        // Partition primitives into two sets and build children
        int mid = (start + end) / 2;
        if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim]) {
            // Create leaf _BVHBuildNode_
            int firstPrimOffset = orderedPrimsOffset->fetch_add(nPrimitives);
            for (int i = start; i < end; ++i) {
                int index = bvhPrimitives[i].primitiveIndex;
                orderedPrims[firstPrimOffset + i - start] = primitives[index];
            }
            node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
            return node;

        } else {
            // Partition primitives based on _splitMethod_
            switch (splitMethod) {
            case SplitMethod::Middle: {
                // Partition primitives through node's midpoint
                Float pmid = (centroidBounds.pMin[dim] + centroidBounds.pMax[dim]) / 2;
                BVHPrimitive *midPtr =
                    std::partition(&bvhPrimitives[start], &bvhPrimitives[end - 1] + 1,
                                   [dim, pmid](const BVHPrimitive &pi) {
                                       return pi.centroid[dim] < pmid;
                                   });
                mid = midPtr - &bvhPrimitives[0];
                // For lots of prims with large overlapping bounding boxes, this
                // may fail to partition; in that case don't break and fall through
                // to EqualCounts.
                if (mid != start && mid != end)
                    break;
            }
            case SplitMethod::EqualCounts: {
                // Partition primitives into equally sized subsets
                mid = (start + end) / 2;
                std::nth_element(&bvhPrimitives[start], &bvhPrimitives[mid],
                                 &bvhPrimitives[end - 1] + 1,
                                 [dim](const BVHPrimitive &a, const BVHPrimitive &b) {
                                     return a.centroid[dim] < b.centroid[dim];
                                 });

                break;
            }
            case SplitMethod::SAH:
            default: {
                // Partition primitives using approximate SAH
                if (nPrimitives <= 2) {
                    // Partition primitives into equally sized subsets
                    mid = (start + end) / 2;
                    std::nth_element(&bvhPrimitives[start], &bvhPrimitives[mid],
                                     &bvhPrimitives[end - 1] + 1,
                                     [dim](const BVHPrimitive &a, const BVHPrimitive &b) {
                                         return a.centroid[dim] < b.centroid[dim];
                                     });

                } else {
                    // Allocate _BVHSplitBucket_ for SAH partition buckets
                    constexpr int nBuckets = 12;
                    BVHSplitBucket buckets[nBuckets];

                    // Initialize _BVHSplitBucket_ for SAH partition buckets
                    for (int i = start; i < end; ++i) {
                        int b = nBuckets *
                                centroidBounds.Offset(bvhPrimitives[i].centroid)[dim];
                        if (b == nBuckets)
                            b = nBuckets - 1;
                        DCHECK_GE(b, 0);
                        DCHECK_LT(b, nBuckets);
                        buckets[b].count++;
                        buckets[b].bounds =
                            Union(buckets[b].bounds, bvhPrimitives[i].bounds);
                    }

                    // Compute costs for splitting after each bucket
                    constexpr int nSplits = nBuckets - 1;
                    int countBelow[nSplits], countAbove[nSplits];
                    Bounds3f boundsBelow[nSplits], boundsAbove[nSplits];
                    // Initialize _countBelow_ and _boundsBelow_ using a forward scan over
                    // splits
                    countBelow[0] = buckets[0].count;
                    boundsBelow[0] = buckets[0].bounds;
                    for (int i = 1; i < nSplits; ++i) {
                        countBelow[i] = countBelow[i - 1] + buckets[i].count;
                        boundsBelow[i] = Union(boundsBelow[i - 1], buckets[i].bounds);
                    }

                    // Initialize _countAbove_ and _boundsAbove_ using a backwards scan
                    // over splits
                    countAbove[nSplits - 1] = buckets[nBuckets - 1].count;
                    boundsAbove[nSplits - 1] = buckets[nBuckets - 1].bounds;
                    for (int i = nSplits - 2; i >= 0; --i) {
                        countAbove[i] = countAbove[i + 1] + buckets[i + 1].count;
                        boundsAbove[i] = Union(boundsAbove[i + 1], buckets[i + 1].bounds);
                    }

                    // Find bucket to split at that minimizes SAH metric
                    int minCostSplitBucket = -1;
                    Float minCost = Infinity;
                    for (int i = 0; i < nSplits; ++i) {
                        // Compute cost for candidate split and update minimum if
                        // necessary
                        if (countBelow[i] == 0 || countAbove[i] == 0)
                            continue;
                        Float cost = (countBelow[i] * boundsBelow[i].SurfaceArea() +
                                      countAbove[i] * boundsAbove[i].SurfaceArea());
                        if (cost < minCost) {
                            minCost = cost;
                            minCostSplitBucket = i;
                        }
                    }
                    // Compute leaf cost and SAH split cost for chosen split
                    Float leafCost = nPrimitives;
                    minCost = 1.f / 2.f + minCost / bounds.SurfaceArea();

                    // Either create leaf or split primitives at selected SAH bucket
                    if (nPrimitives > maxPrimsInNode || minCost < leafCost) {
                        BVHPrimitive *pmid = std::partition(
                            &bvhPrimitives[start], &bvhPrimitives[end - 1] + 1,
                            [=](const BVHPrimitive &bp) {
                                int b =
                                    nBuckets * centroidBounds.Offset(bp.centroid)[dim];
                                if (b == nBuckets)
                                    b = nBuckets - 1;
                                return b <= minCostSplitBucket;
                            });
                        mid = pmid - &bvhPrimitives[0];
                    } else {
                        // Create leaf _BVHBuildNode_
                        int firstPrimOffset = orderedPrimsOffset->fetch_add(nPrimitives);
                        for (int i = start; i < end; ++i) {
                            int index = bvhPrimitives[i].primitiveIndex;
                            orderedPrims[firstPrimOffset + i - start] = primitives[index];
                        }
                        node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
                        return node;
                    }
                }

                break;
            }
            }

            BVHBuildNode *children[2];
            // Recursively build BVHs for _children_
            if (end - start > 128 * 1024) {
                // Recursively build child BVHs in parallel
                ParallelFor(0, 2, [&](int i) {
                    if (i == 0)
                        children[0] =
                            buildRecursive(threadAllocators, bvhPrimitives, start, mid,
                                           totalNodes, orderedPrimsOffset, orderedPrims);
                    else
                        children[1] =
                            buildRecursive(threadAllocators, bvhPrimitives, mid, end,
                                           totalNodes, orderedPrimsOffset, orderedPrims);
                });

            } else {
                // Recursively build child BVHs sequentially
                children[0] =
                    buildRecursive(threadAllocators, bvhPrimitives, start, mid,
                                   totalNodes, orderedPrimsOffset, orderedPrims);
                children[1] =
                    buildRecursive(threadAllocators, bvhPrimitives, mid, end, totalNodes,
                                   orderedPrimsOffset, orderedPrims);
            }

            node->InitInterior(dim, children[0], children[1]);
        }
    }

    return node;
}

BVHBuildNode *BVHAggregate::buildHLBVH(Allocator alloc,
                                       const std::vector<BVHPrimitive> &bvhPrimitives,
                                       std::atomic<int> *totalNodes,
                                       std::vector<PrimitiveHandle> &orderedPrims) {
    // Compute bounding box of all primitive centroids
    Bounds3f bounds;
    for (const BVHPrimitive &pi : bvhPrimitives)
        bounds = Union(bounds, pi.centroid);

    // Compute Morton indices of primitives
    std::vector<MortonPrimitive> mortonPrims(bvhPrimitives.size());
    ParallelFor(0, bvhPrimitives.size(), [&](int64_t i) {
        // Initialize _mortonPrims[i]_ for _i_th primitive
        constexpr int mortonBits = 10;
        constexpr int mortonScale = 1 << mortonBits;
        mortonPrims[i].primitiveIndex = bvhPrimitives[i].primitiveIndex;
        Vector3f centroidOffset = bounds.Offset(bvhPrimitives[i].centroid);
        Vector3f offset = centroidOffset * mortonScale;
        mortonPrims[i].mortonCode = EncodeMorton3(offset.x, offset.y, offset.z);
    });

    // Radix sort primitive Morton indices
    RadixSort(&mortonPrims);

    // Create LBVH treelets at bottom of BVH
    // Find intervals of primitives for each treelet
    std::vector<LBVHTreelet> treeletsToBuild;
    for (int start = 0, end = 1; end <= (int)mortonPrims.size(); ++end) {
        uint32_t mask = 0b00111111111111000000000000000000;
        if (end == (int)mortonPrims.size() || ((mortonPrims[start].mortonCode & mask) !=
                                               (mortonPrims[end].mortonCode & mask))) {
            // Add entry to _treeletsToBuild_ for this treelet
            int nPrimitives = end - start;
            int maxBVHNodes = 2 * nPrimitives - 1;
            BVHBuildNode *nodes = alloc.allocate_object<BVHBuildNode>(maxBVHNodes);
            treeletsToBuild.push_back({start, nPrimitives, nodes});

            start = end;
        }
    }

    // Create LBVHs for treelets in parallel
    std::atomic<int> orderedPrimsOffset(0);
    ParallelFor(0, treeletsToBuild.size(), [&](int i) {
        // Generate _i_th LBVH treelet
        int nodesCreated = 0;
        const int firstBitIndex = 29 - 12;
        LBVHTreelet &tr = treeletsToBuild[i];
        tr.buildNodes = emitLBVH(
            tr.buildNodes, bvhPrimitives, &mortonPrims[tr.startIndex], tr.nPrimitives,
            &nodesCreated, orderedPrims, &orderedPrimsOffset, firstBitIndex);
        *totalNodes += nodesCreated;
    });

    // Create and return SAH BVH from LBVH treelets
    std::vector<BVHBuildNode *> finishedTreelets;
    finishedTreelets.reserve(treeletsToBuild.size());
    for (LBVHTreelet &treelet : treeletsToBuild)
        finishedTreelets.push_back(treelet.buildNodes);
    return buildUpperSAH(alloc, finishedTreelets, 0, finishedTreelets.size(), totalNodes);
}

BVHBuildNode *BVHAggregate::emitLBVH(BVHBuildNode *&buildNodes,
                                     const std::vector<BVHPrimitive> &bvhPrimitives,
                                     MortonPrimitive *mortonPrims, int nPrimitives,
                                     int *totalNodes,
                                     std::vector<PrimitiveHandle> &orderedPrims,
                                     std::atomic<int> *orderedPrimsOffset, int bitIndex) {
    CHECK_GT(nPrimitives, 0);
    if (bitIndex == -1 || nPrimitives < maxPrimsInNode) {
        // Create and return leaf node of LBVH treelet
        ++*totalNodes;
        BVHBuildNode *node = buildNodes++;
        Bounds3f bounds;
        int firstPrimOffset = orderedPrimsOffset->fetch_add(nPrimitives);
        for (int i = 0; i < nPrimitives; ++i) {
            int primitiveIndex = mortonPrims[i].primitiveIndex;
            orderedPrims[firstPrimOffset + i] = primitives[primitiveIndex];
            bounds = Union(bounds, bvhPrimitives[primitiveIndex].bounds);
        }
        node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
        return node;

    } else {
        int mask = 1 << bitIndex;
        // Advance to next subtree level if there's no LBVH split for this bit
        if ((mortonPrims[0].mortonCode & mask) ==
            (mortonPrims[nPrimitives - 1].mortonCode & mask))
            return emitLBVH(buildNodes, bvhPrimitives, mortonPrims, nPrimitives,
                            totalNodes, orderedPrims, orderedPrimsOffset, bitIndex - 1);

        // Find LBVH split point for this dimension
        int splitOffset = FindInterval(nPrimitives, [&](int index) {
            return ((mortonPrims[0].mortonCode & mask) ==
                    (mortonPrims[index].mortonCode & mask));
        });
        ++splitOffset;
        CHECK_LE(splitOffset, nPrimitives - 1);
        CHECK_NE(mortonPrims[splitOffset - 1].mortonCode & mask,
                 mortonPrims[splitOffset].mortonCode & mask);

        // Create and return interior LBVH node
        (*totalNodes)++;
        BVHBuildNode *node = buildNodes++;
        BVHBuildNode *lbvh[2] = {
            emitLBVH(buildNodes, bvhPrimitives, mortonPrims, splitOffset, totalNodes,
                     orderedPrims, orderedPrimsOffset, bitIndex - 1),
            emitLBVH(buildNodes, bvhPrimitives, &mortonPrims[splitOffset],
                     nPrimitives - splitOffset, totalNodes, orderedPrims,
                     orderedPrimsOffset, bitIndex - 1)};
        int axis = bitIndex % 3;
        node->InitInterior(axis, lbvh[0], lbvh[1]);
        return node;
    }
}

int BVHAggregate::flattenBVHTree(BVHBuildNode *node, int *offset) {
    LinearBVHNode *linearNode = &nodes[*offset];
    linearNode->bounds = node->bounds;
    int nodeOffset = (*offset)++;
    if (node->nPrimitives > 0) {
        CHECK(!node->children[0] && !node->children[1]);
        CHECK_LT(node->nPrimitives, 65536);
        linearNode->primitivesOffset = node->firstPrimOffset;
        linearNode->nPrimitives = node->nPrimitives;
    } else {
        // Create interior flattened BVH node
        linearNode->axis = node->splitAxis;
        linearNode->nPrimitives = 0;
        flattenBVHTree(node->children[0], offset);
        linearNode->secondChildOffset = flattenBVHTree(node->children[1], offset);
    }
    return nodeOffset;
}

Bounds3f BVHAggregate::Bounds() const {
    CHECK(nodes != nullptr);
    return nodes[0].bounds;
}

pstd::optional<ShapeIntersection> BVHAggregate::Intersect(const Ray &ray,
                                                          Float tMax) const {
    if (nodes == nullptr)
        return {};
    pstd::optional<ShapeIntersection> si;
    Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
    int dirIsNeg[3] = {int(invDir.x < 0), int(invDir.y < 0), int(invDir.z < 0)};
    // Follow ray through BVH nodes to find primitive intersections
    int toVisitOffset = 0, currentNodeIndex = 0;
    int nodesToVisit[64];
    int nodesVisited = 0;
    while (true) {
        ++nodesVisited;
        const LinearBVHNode *node = &nodes[currentNodeIndex];
        // Check ray against BVH node
        if (node->bounds.IntersectP(ray.o, ray.d, tMax, invDir, dirIsNeg)) {
            if (node->nPrimitives > 0) {
                // Intersect ray with primitives in leaf BVH node
                for (int i = 0; i < node->nPrimitives; ++i) {
                    // Check for intersection with primitive in BVH node
                    pstd::optional<ShapeIntersection> primSi =
                        primitives[node->primitivesOffset + i].Intersect(ray, tMax);
                    if (primSi) {
                        si = primSi;
                        tMax = si->tHit;
                    }
                }
                if (toVisitOffset == 0)
                    break;
                currentNodeIndex = nodesToVisit[--toVisitOffset];

            } else {
                // Put far BVH node on _nodesToVisit_ stack, advance to near node
                if (dirIsNeg[node->axis]) {
                    nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
                    currentNodeIndex = node->secondChildOffset;
                } else {
                    nodesToVisit[toVisitOffset++] = node->secondChildOffset;
                    currentNodeIndex = currentNodeIndex + 1;
                }
            }
        } else {
            if (toVisitOffset == 0)
                break;
            currentNodeIndex = nodesToVisit[--toVisitOffset];
        }
    }

    bvhNodesVisited += nodesVisited;
    return si;
}

bool BVHAggregate::IntersectP(const Ray &ray, Float tMax) const {
    if (nodes == nullptr)
        return false;
    Vector3f invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
    int dirIsNeg[3] = {static_cast<int>(invDir.x < 0), static_cast<int>(invDir.y < 0),
                       static_cast<int>(invDir.z < 0)};
    int nodesToVisit[64];
    int toVisitOffset = 0, currentNodeIndex = 0;
    int nodesVisited = 0;

    while (true) {
        ++nodesVisited;
        const LinearBVHNode *node = &nodes[currentNodeIndex];
        if (node->bounds.IntersectP(ray.o, ray.d, tMax, invDir, dirIsNeg)) {
            // Process BVH node _node_ for traversal
            if (node->nPrimitives > 0) {
                for (int i = 0; i < node->nPrimitives; ++i) {
                    if (primitives[node->primitivesOffset + i].IntersectP(ray, tMax)) {
                        bvhNodesVisited += nodesVisited;
                        return true;
                    }
                }
                if (toVisitOffset == 0)
                    break;
                currentNodeIndex = nodesToVisit[--toVisitOffset];
            } else {
                if (dirIsNeg[node->axis] != 0) {
                    /// second child first
                    nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
                    currentNodeIndex = node->secondChildOffset;
                } else {
                    nodesToVisit[toVisitOffset++] = node->secondChildOffset;
                    currentNodeIndex = currentNodeIndex + 1;
                }
            }
        } else {
            if (toVisitOffset == 0)
                break;
            currentNodeIndex = nodesToVisit[--toVisitOffset];
        }
    }
    bvhNodesVisited += nodesVisited;
    return false;
}

BVHBuildNode *BVHAggregate::buildUpperSAH(Allocator alloc,
                                          std::vector<BVHBuildNode *> &treeletRoots,
                                          int start, int end,
                                          std::atomic<int> *totalNodes) const {
    CHECK_LT(start, end);
    int nNodes = end - start;
    if (nNodes == 1)
        return treeletRoots[start];
    (*totalNodes)++;
    BVHBuildNode *node = alloc.new_object<BVHBuildNode>();

    // Compute bounds of all nodes under this HLBVH node
    Bounds3f bounds;
    for (int i = start; i < end; ++i)
        bounds = Union(bounds, treeletRoots[i]->bounds);

    // Compute bound of HLBVH node centroids, choose split dimension _dim_
    Bounds3f centroidBounds;
    for (int i = start; i < end; ++i) {
        Point3f centroid =
            (treeletRoots[i]->bounds.pMin + treeletRoots[i]->bounds.pMax) * 0.5f;
        centroidBounds = Union(centroidBounds, centroid);
    }
    int dim = centroidBounds.MaxDimension();
    // FIXME: if this hits, what do we need to do?
    // Make sure the SAH split below does something... ?
    CHECK_NE(centroidBounds.pMax[dim], centroidBounds.pMin[dim]);

    // Allocate _BVHSplitBucket_ for SAH partition buckets
    constexpr int nBuckets = 12;
    struct BVHSplitBucket {
        int count = 0;
        Bounds3f bounds;
    };
    BVHSplitBucket buckets[nBuckets];

    // Initialize _BVHSplitBucket_ for HLBVH SAH partition buckets
    for (int i = start; i < end; ++i) {
        Float centroid =
            (treeletRoots[i]->bounds.pMin[dim] + treeletRoots[i]->bounds.pMax[dim]) *
            0.5f;
        int b = nBuckets * ((centroid - centroidBounds.pMin[dim]) /
                            (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
        if (b == nBuckets)
            b = nBuckets - 1;
        CHECK_GE(b, 0);
        CHECK_LT(b, nBuckets);
        buckets[b].count++;
        buckets[b].bounds = Union(buckets[b].bounds, treeletRoots[i]->bounds);
    }

    // Compute costs for splitting after each bucket
    Float cost[nBuckets - 1];
    for (int i = 0; i < nBuckets - 1; ++i) {
        Bounds3f b0, b1;
        int count0 = 0, count1 = 0;
        for (int j = 0; j <= i; ++j) {
            b0 = Union(b0, buckets[j].bounds);
            count0 += buckets[j].count;
        }
        for (int j = i + 1; j < nBuckets; ++j) {
            b1 = Union(b1, buckets[j].bounds);
            count1 += buckets[j].count;
        }
        cost[i] = .125f + (count0 * b0.SurfaceArea() + count1 * b1.SurfaceArea()) /
                              bounds.SurfaceArea();
    }

    // Find bucket to split at that minimizes SAH metric
    Float minCost = cost[0];
    int minCostSplitBucket = 0;
    for (int i = 1; i < nBuckets - 1; ++i) {
        if (cost[i] < minCost) {
            minCost = cost[i];
            minCostSplitBucket = i;
        }
    }

    // Split nodes and create interior HLBVH SAH node
    BVHBuildNode **pmid = std::partition(
        &treeletRoots[start], &treeletRoots[end - 1] + 1, [=](const BVHBuildNode *node) {
            Float centroid = (node->bounds.pMin[dim] + node->bounds.pMax[dim]) * 0.5f;
            int b = nBuckets * ((centroid - centroidBounds.pMin[dim]) /
                                (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
            if (b == nBuckets)
                b = nBuckets - 1;
            CHECK_GE(b, 0);
            CHECK_LT(b, nBuckets);
            return b <= minCostSplitBucket;
        });
    int mid = pmid - &treeletRoots[0];
    CHECK_GT(mid, start);
    CHECK_LT(mid, end);
    node->InitInterior(dim,
                       this->buildUpperSAH(alloc, treeletRoots, start, mid, totalNodes),
                       this->buildUpperSAH(alloc, treeletRoots, mid, end, totalNodes));
    return node;
}

BVHAggregate *BVHAggregate::Create(std::vector<PrimitiveHandle> prims,
                                   const ParameterDictionary &parameters) {
    std::string splitMethodName = parameters.GetOneString("splitmethod", "sah");
    BVHAggregate::SplitMethod splitMethod;
    if (splitMethodName == "sah")
        splitMethod = BVHAggregate::SplitMethod::SAH;
    else if (splitMethodName == "hlbvh")
        splitMethod = BVHAggregate::SplitMethod::HLBVH;
    else if (splitMethodName == "middle")
        splitMethod = BVHAggregate::SplitMethod::Middle;
    else if (splitMethodName == "equal")
        splitMethod = BVHAggregate::SplitMethod::EqualCounts;
    else {
        Warning(R"(BVH split method "%s" unknown.  Using "sah".)", splitMethodName);
        splitMethod = BVHAggregate::SplitMethod::SAH;
    }

    int maxPrimsInNode = parameters.GetOneInt("maxnodeprims", 4);
    return new BVHAggregate(std::move(prims), maxPrimsInNode, splitMethod);
}

// KdNodeToVisit Definition
struct KdNodeToVisit {
    const KdTreeNode *node;
    Float tMin, tMax;
};

// KdTreeNode Definition
struct alignas(8) KdTreeNode {
    // KdTreeNode Methods
    void InitLeaf(int *primNums, int np, std::vector<int> *primitiveIndices);

    void InitInterior(int axis, int ac, Float s) {
        split = s;
        flags = axis;
        aboveChild |= (ac << 2);
    }

    Float SplitPos() const { return split; }
    int nPrimitives() const { return nPrims >> 2; }
    int SplitAxis() const { return flags & 3; }
    bool IsLeaf() const { return (flags & 3) == 3; }
    int AboveChild() const { return aboveChild >> 2; }

    union {
        Float split;                 // Interior
        int onePrimitive;            // Leaf
        int primitiveIndicesOffset;  // Leaf
    };

  private:
    union {
        int flags;       // Both
        int nPrims;      // Leaf
        int aboveChild;  // Interior
    };
};

// EdgeType Definition
enum class EdgeType { Start, End };

// BoundEdge Definition
struct BoundEdge {
    // BoundEdge Public Methods
    BoundEdge() {}

    BoundEdge(Float t, int primNum, bool starting) : t(t), primNum(primNum) {
        type = starting ? EdgeType::Start : EdgeType::End;
    }

    Float t;
    int primNum;
    EdgeType type;
};

STAT_PIXEL_COUNTER("Kd-Tree/Nodes visited", kdNodesVisited);

// KdTreeAggregate Method Definitions
KdTreeAggregate::KdTreeAggregate(std::vector<PrimitiveHandle> p, int isectCost,
                                 int traversalCost, Float emptyBonus, int maxPrims,
                                 int maxDepth)
    : isectCost(isectCost),
      traversalCost(traversalCost),
      maxPrims(maxPrims),
      emptyBonus(emptyBonus),
      primitives(std::move(p)) {
    // Build kd-tree for accelerator
    nextFreeNode = nAllocedNodes = 0;
    if (maxDepth <= 0)
        maxDepth = std::round(8 + 1.3f * Log2Int(int64_t(primitives.size())));
    // Compute bounds for kd-tree construction
    std::vector<Bounds3f> primBounds;
    primBounds.reserve(primitives.size());
    for (PrimitiveHandle &prim : primitives) {
        Bounds3f b = prim.Bounds();
        bounds = Union(bounds, b);
        primBounds.push_back(b);
    }

    // Allocate working memory for kd-tree construction
    std::unique_ptr<BoundEdge[]> edges[3];
    for (int i = 0; i < 3; ++i)
        edges[i] = std::make_unique<BoundEdge[]>(2 * primitives.size());

    std::unique_ptr<int[]> prims0 = std::make_unique<int[]>(primitives.size());
    std::unique_ptr<int[]> prims1 =
        std::make_unique<int[]>((maxDepth + 1) * primitives.size());

    // Initialize _primNums_ for kd-tree construction
    std::unique_ptr<int[]> primNums = std::make_unique<int[]>(primitives.size());
    for (size_t i = 0; i < primitives.size(); ++i)
        primNums[i] = i;

    // Start recursive construction of kd-tree
    buildTree(0, bounds, primBounds, primNums.get(), primitives.size(), maxDepth, edges,
              prims0.get(), prims1.get(), 0);
}

void KdTreeNode::InitLeaf(int *primNums, int np, std::vector<int> *primitiveIndices) {
    flags = 3;
    nPrims |= (np << 2);
    // Store primitive ids for leaf node
    if (np == 0)
        onePrimitive = 0;
    else if (np == 1)
        onePrimitive = primNums[0];
    else {
        primitiveIndicesOffset = primitiveIndices->size();
        for (int i = 0; i < np; ++i)
            primitiveIndices->push_back(primNums[i]);
    }
}

void KdTreeAggregate::buildTree(int nodeNum, const Bounds3f &nodeBounds,
                                const std::vector<Bounds3f> &allPrimBounds, int *primNums,
                                int nPrimitives, int depth,
                                const std::unique_ptr<BoundEdge[]> edges[3], int *prims0,
                                int *prims1, int badRefines) {
    CHECK_EQ(nodeNum, nextFreeNode);
    // Get next free node from _nodes_ array
    if (nextFreeNode == nAllocedNodes) {
        int nNewAllocNodes = std::max(2 * nAllocedNodes, 512);
        KdTreeNode *n = new KdTreeNode[nNewAllocNodes];
        if (nAllocedNodes > 0) {
            memcpy(n, nodes, nAllocedNodes * sizeof(KdTreeNode));
            delete[] nodes;
        }
        nodes = n;
        nAllocedNodes = nNewAllocNodes;
    }
    ++nextFreeNode;

    // Initialize leaf node if termination criteria met
    if (nPrimitives <= maxPrims || depth == 0) {
        nodes[nodeNum].InitLeaf(primNums, nPrimitives, &primitiveIndices);
        return;
    }

    // Initialize interior node and continue recursion
    // Choose split axis position for interior node
    int bestAxis = -1, bestOffset = -1;
    Float bestCost = Infinity, leafCost = isectCost * nPrimitives;
    Float invTotalSA = 1 / nodeBounds.SurfaceArea();
    // Choose which axis to split along
    int axis = nodeBounds.MaxDimension();

    // Choose split along axis and attempt to partition primitives
    int retries = 0;
retrySplit:
    // Initialize edges for _axis_
    for (int i = 0; i < nPrimitives; ++i) {
        int pn = primNums[i];
        const Bounds3f &bounds = allPrimBounds[pn];
        edges[axis][2 * i] = BoundEdge(bounds.pMin[axis], pn, true);
        edges[axis][2 * i + 1] = BoundEdge(bounds.pMax[axis], pn, false);
    }
    // Sort _edges_ for _axis_
    std::sort(&edges[axis][0], &edges[axis][2 * nPrimitives],
              [](const BoundEdge &e0, const BoundEdge &e1) -> bool {
                  return std::tie(e0.t, e0.type) < std::tie(e1.t, e1.type);
              });

    // Compute cost of all splits for _axis_ to find best
    int nBelow = 0, nAbove = nPrimitives;
    for (int i = 0; i < 2 * nPrimitives; ++i) {
        if (edges[axis][i].type == EdgeType::End)
            --nAbove;
        Float edgeT = edges[axis][i].t;
        if (edgeT > nodeBounds.pMin[axis] && edgeT < nodeBounds.pMax[axis]) {
            // Compute child surface areas for split at _edgeT_
            Vector3f d = nodeBounds.pMax - nodeBounds.pMin;
            int otherAxis0 = (axis + 1) % 3, otherAxis1 = (axis + 2) % 3;
            Float belowSA =
                2 * (d[otherAxis0] * d[otherAxis1] +
                     (edgeT - nodeBounds.pMin[axis]) * (d[otherAxis0] + d[otherAxis1]));
            Float aboveSA =
                2 * (d[otherAxis0] * d[otherAxis1] +
                     (nodeBounds.pMax[axis] - edgeT) * (d[otherAxis0] + d[otherAxis1]));

            // Compute cost for split at _i_th edge
            Float pBelow = belowSA * invTotalSA, pAbove = aboveSA * invTotalSA;
            Float eb = (nAbove == 0 || nBelow == 0) ? emptyBonus : 0;
            Float cost = traversalCost +
                         isectCost * (1 - eb) * (pBelow * nBelow + pAbove * nAbove);
            // Update best split if this is lowest cost so far
            if (cost < bestCost) {
                bestCost = cost;
                bestAxis = axis;
                bestOffset = i;
            }
        }
        if (edges[axis][i].type == EdgeType::Start)
            ++nBelow;
    }
    CHECK(nBelow == nPrimitives && nAbove == 0);

    // Create leaf if no good splits were found
    if (bestAxis == -1 && retries < 2) {
        ++retries;
        axis = (axis + 1) % 3;
        goto retrySplit;
    }
    if (bestCost > leafCost)
        ++badRefines;
    if ((bestCost > 4 * leafCost && nPrimitives < 16) || bestAxis == -1 ||
        badRefines == 3) {
        nodes[nodeNum].InitLeaf(primNums, nPrimitives, &primitiveIndices);
        return;
    }

    // Classify primitives with respect to split
    int n0 = 0, n1 = 0;
    for (int i = 0; i < bestOffset; ++i)
        if (edges[bestAxis][i].type == EdgeType::Start)
            prims0[n0++] = edges[bestAxis][i].primNum;
    for (int i = bestOffset + 1; i < 2 * nPrimitives; ++i)
        if (edges[bestAxis][i].type == EdgeType::End)
            prims1[n1++] = edges[bestAxis][i].primNum;

    // Recursively initialize children nodes
    Float tSplit = edges[bestAxis][bestOffset].t;
    Bounds3f bounds0 = nodeBounds, bounds1 = nodeBounds;
    bounds0.pMax[bestAxis] = bounds1.pMin[bestAxis] = tSplit;
    buildTree(nodeNum + 1, bounds0, allPrimBounds, prims0, n0, depth - 1, edges, prims0,
              prims1 + nPrimitives, badRefines);
    int aboveChild = nextFreeNode;
    nodes[nodeNum].InitInterior(bestAxis, aboveChild, tSplit);
    buildTree(aboveChild, bounds1, allPrimBounds, prims1, n1, depth - 1, edges, prims0,
              prims1 + nPrimitives, badRefines);
}

pstd::optional<ShapeIntersection> KdTreeAggregate::Intersect(const Ray &ray,
                                                             Float rayTMax) const {
    // Compute initial parametric range of ray inside kd-tree extent
    Float tMin, tMax;
    if (!bounds.IntersectP(ray.o, ray.d, rayTMax, &tMin, &tMax))
        return {};

    // Prepare to traverse kd-tree for ray
    Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
    constexpr int maxToVisit = 64;
    KdNodeToVisit toVisit[maxToVisit];
    int toVisitIndex = 0;
    int nodesVisited = 0;

    // Traverse kd-tree nodes in order for ray
    pstd::optional<ShapeIntersection> si;
    const KdTreeNode *node = &nodes[0];
    while (node != nullptr) {
        // Bail out if we found a hit closer than the current node
        if (rayTMax < tMin)
            break;

        ++nodesVisited;
        if (!node->IsLeaf()) {
            // Visit kd-tree interior node
            // Compute parametric distance along ray to split plane
            int axis = node->SplitAxis();
            Float tSplit = (node->SplitPos() - ray.o[axis]) * invDir[axis];

            // Get node children pointers for ray
            const KdTreeNode *firstChild, *secondChild;
            int belowFirst = (ray.o[axis] < node->SplitPos()) ||
                             (ray.o[axis] == node->SplitPos() && ray.d[axis] <= 0);
            if (belowFirst) {
                firstChild = node + 1;
                secondChild = &nodes[node->AboveChild()];
            } else {
                firstChild = &nodes[node->AboveChild()];
                secondChild = node + 1;
            }

            // Advance to next child node, possibly enqueue other child
            if (tSplit > tMax || tSplit <= 0)
                node = firstChild;
            else if (tSplit < tMin)
                node = secondChild;
            else {
                // Enqueue _secondChild_ in todo list
                toVisit[toVisitIndex].node = secondChild;
                toVisit[toVisitIndex].tMin = tSplit;
                toVisit[toVisitIndex].tMax = tMax;
                ++toVisitIndex;

                node = firstChild;
                tMax = tSplit;
            }

        } else {
            // Check for intersections inside leaf node
            int nPrimitives = node->nPrimitives();
            if (nPrimitives == 1) {
                const PrimitiveHandle &p = primitives[node->onePrimitive];
                // Check one primitive inside leaf node
                pstd::optional<ShapeIntersection> primSi = p.Intersect(ray, rayTMax);
                if (primSi) {
                    si = primSi;
                    rayTMax = si->tHit;
                }

            } else {
                for (int i = 0; i < nPrimitives; ++i) {
                    int index = primitiveIndices[node->primitiveIndicesOffset + i];
                    const PrimitiveHandle &p = primitives[index];
                    // Check one primitive inside leaf node
                    pstd::optional<ShapeIntersection> primSi = p.Intersect(ray, rayTMax);
                    if (primSi) {
                        si = primSi;
                        rayTMax = si->tHit;
                    }
                }
            }

            // Grab next node to visit from todo list
            if (toVisitIndex > 0) {
                --toVisitIndex;
                node = toVisit[toVisitIndex].node;
                tMin = toVisit[toVisitIndex].tMin;
                tMax = toVisit[toVisitIndex].tMax;
            } else
                break;
        }
    }
    kdNodesVisited += nodesVisited;
    return si;
}

bool KdTreeAggregate::IntersectP(const Ray &ray, Float raytMax) const {
    // Compute initial parametric range of ray inside kd-tree extent
    Float tMin, tMax;
    if (!bounds.IntersectP(ray.o, ray.d, raytMax, &tMin, &tMax))
        return false;

    // Prepare to traverse kd-tree for ray
    Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
    constexpr int maxTodo = 64;
    KdNodeToVisit toVisit[maxTodo];
    int toVisitIndex = 0;
    int nodesVisited = 0;
    const KdTreeNode *node = &nodes[0];
    while (node != nullptr) {
        ++nodesVisited;
        if (node->IsLeaf()) {
            // Check for shadow ray intersections inside leaf node
            int nPrimitives = node->nPrimitives();
            if (nPrimitives == 1) {
                const PrimitiveHandle &p = primitives[node->onePrimitive];
                if (p.IntersectP(ray, raytMax)) {
                    kdNodesVisited += nodesVisited;
                    return true;
                }
            } else {
                for (int i = 0; i < nPrimitives; ++i) {
                    int primitiveIndex =
                        primitiveIndices[node->primitiveIndicesOffset + i];
                    const PrimitiveHandle &prim = primitives[primitiveIndex];
                    if (prim.IntersectP(ray, raytMax)) {
                        kdNodesVisited += nodesVisited;
                        return true;
                    }
                }
            }

            // Grab next node to process from todo list
            if (toVisitIndex > 0) {
                --toVisitIndex;
                node = toVisit[toVisitIndex].node;
                tMin = toVisit[toVisitIndex].tMin;
                tMax = toVisit[toVisitIndex].tMax;
            } else
                break;
        } else {
            // Process kd-tree interior node

            // Compute parametric distance along ray to split plane
            int axis = node->SplitAxis();
            Float tSplit = (node->SplitPos() - ray.o[axis]) * invDir[axis];

            // Get node children pointers for ray
            const KdTreeNode *firstChild, *secondChild;
            int belowFirst = (ray.o[axis] < node->SplitPos()) ||
                             (ray.o[axis] == node->SplitPos() && ray.d[axis] <= 0);
            if (belowFirst != 0) {
                firstChild = node + 1;
                secondChild = &nodes[node->AboveChild()];
            } else {
                firstChild = &nodes[node->AboveChild()];
                secondChild = node + 1;
            }

            // Advance to next child node, possibly enqueue other child
            if (tSplit > tMax || tSplit <= 0)
                node = firstChild;
            else if (tSplit < tMin)
                node = secondChild;
            else {
                // Enqueue _secondChild_ in todo list
                toVisit[toVisitIndex].node = secondChild;
                toVisit[toVisitIndex].tMin = tSplit;
                toVisit[toVisitIndex].tMax = tMax;
                ++toVisitIndex;
                node = firstChild;
                tMax = tSplit;
            }
        }
    }
    kdNodesVisited += nodesVisited;
    return false;
}

KdTreeAggregate *KdTreeAggregate::Create(std::vector<PrimitiveHandle> prims,
                                         const ParameterDictionary &parameters) {
    int isectCost = parameters.GetOneInt("intersectcost", 5);
    int travCost = parameters.GetOneInt("traversalcost", 1);
    Float emptyBonus = parameters.GetOneFloat("emptybonus", 0.5f);
    int maxPrims = parameters.GetOneInt("maxprims", 1);
    int maxDepth = parameters.GetOneInt("maxdepth", -1);
    return new KdTreeAggregate(std::move(prims), isectCost, travCost, emptyBonus,
                               maxPrims, maxDepth);
}

PrimitiveHandle CreateAccelerator(const std::string &name,
                                  std::vector<PrimitiveHandle> prims,
                                  const ParameterDictionary &parameters) {
    PrimitiveHandle accel = nullptr;
    if (name == "bvh")
        accel = BVHAggregate::Create(std::move(prims), parameters);
    else if (name == "kdtree")
        accel = KdTreeAggregate::Create(std::move(prims), parameters);
    else
        ErrorExit("%s: accelerator type unknown.", name);

    if (!accel)
        ErrorExit("%s: unable to create accelerator.", name);

    parameters.ReportUnused();
    return accel;
}

}  // namespace pbrt
