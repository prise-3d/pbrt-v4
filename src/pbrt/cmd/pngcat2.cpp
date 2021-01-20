#include <cstdio>
#include <png.h>
#include <iostream>
using namespace std;



  static const int LOADED=0; /**<Constante permettant de signifier que le chargement d'une image s'est déroulé sans aucun problème.*/
  static const int ERR_OPENING=-1; /**<Constante permettant de signifier qu'une erreur s'est produite lors de l'ouverture du fichier image spécifié.*/
  static const int ERR_READING=-2; /**<Constante permettant de signifier qu'une erreur s'est produite lors de la lecture du fichier image spécifié.*/
  static const int ERR_BADFORMAT=-3; /**<Constante permettant de signifier que le format de l'image contenu dans le fichier spécifié n'est pas reconnu par le "loader".*/
  static const int ERR_MEMORY=-4; /**<Constante permettant de signifier que l'image n'a pas pu être allouée en mémoire.*/

int writePng(char *file_name, png_byte *img,
	     int width, int height, int channel);

int readPng(char *filename, int *width,
	    int *height, int *channels, unsigned char **data);



int main(int argc , char *argv[]){
  
  // test des paramètres
  if(argc!=4){
    cout << "syntaxe : " << argv[0] << " fileleft.png fileright.png  file.png" << endl;
    return 0;
  }

  // chargement des deux images
  int wl, hl, wr, hr;// dimensions des deux images left et right
  int cl, cr; // nb de canaux de couleurs pour chaque image
  unsigned char *dl=0, *dr=0; // zones de données pour chaque image
 

  if(readPng(argv[1], &wl, &hl, &cl, &dl)!=LOADED){
    cout << "error while loading " << argv[1] << endl;
    return 0;
  }

  if(readPng(argv[2], &wr, &hr, &cr, &dr)!=LOADED){
    cout << "error while loading " << argv[2] << endl;
    return 0;
  }



  // vérification de la compatibilité des dimensions des deux images
  // On impose que les images aient les mêmes dimensions
  if((wl!=wr) || (hl!=hr) || (cl!=cr)){
    cout << "incompatible image sizes : (";
    cout << wl << "x" << hl << ")-(";
    cout << wr << "x" << hr << "),(";
    cout << cl << "," << cr << ") channels" << endl;
    delete [] dl; delete [] dr;
    return 0;
  }

  cout << "image sizes : (";
  cout << wl << "x" << hl << ")-(";
  cout << wr << "x" << hr << "),(";
  cout << cl << "," << cr << ") channels" << endl;

  // création de l'image de largeur double aux
  // images à y insérer
  int ws = wl*2, hs = hl, cs = cl;
  unsigned char* ds = new unsigned char[ws*hs*cs];

  int il=0, ir=0;

  for(int y=0; y<hs; y++){
    int i=((hs-1) -y)*ws*cs ;
    for(int x=0; x<ws; x++){
      for(int c=0; c<cs; c++){
	ds[i++] = (x<wl) ? dl[il++] : dr[ir++];
      }
    }
  }
  
  delete [] dl; delete [] dr;
  

  // sauvegarde de l'image doublée
  cout << "save" << endl;
  writePng(argv[3], (png_byte*)ds, ws, hs, cs);

  delete [] ds;

  return 1;
}


/* write a png file */
int writePng(char *file_name, png_byte *img, int width, int height, int channel)
{
  FILE *fp;
  png_structp png_ptr;
  png_infop info_ptr;
  //   png_colorp palette;

  cout << " sauvegarde : " << width << "x" << height << "x" << channel << endl;
  /* open the file */
  fp = fopen(file_name, "wb");
  if (fp == NULL)
    return (1);

  /* Create and initialize the png_struct with the desired error handler
   * functions.  If you want to use the default stderr and longjump method,
   * you can supply NULL for the last three parameters.  We also check that
   * the library version is compatible with the one used at compile time,
   * in case we are using dynamically linked libraries.  REQUIRED.
   */
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,NULL, NULL, NULL);

  if (png_ptr == NULL)
    {
      fclose(fp);
      return (1);
    }

  /* Allocate/initialize the image information data.  REQUIRED */
  info_ptr = png_create_info_struct(png_ptr);
  if (info_ptr == NULL)
    {
      fclose(fp);
      png_destroy_write_struct(&png_ptr, 0);
      return (1);
    }

  /* Set error handling.  REQUIRED if you aren't supplying your own
   * error handling functions in the png_create_write_struct() call.
   */
  if (setjmp(png_jmpbuf(png_ptr)))
    {
      /* If we get here, we had a problem reading the file */
      fclose(fp);
      png_destroy_write_struct(&png_ptr, &info_ptr);
      return (1);
    }

  /* set up the output control if you are using standard C streams */
  png_init_io(png_ptr, fp);

  /* This is the hard way */

  /* Set the image information here.  Width and height are up to 2^31,
   * bit_depth is one of 1, 2, 4, 8, or 16, but valid values also depend on
   * the color_type selected. color_type is one of PNG_COLOR_TYPE_GRAY,
   * PNG_COLOR_TYPE_GRAY_ALPHA, PNG_COLOR_TYPE_PALETTE, PNG_COLOR_TYPE_RGB,
   * or PNG_COLOR_TYPE_RGB_ALPHA.  interlace is either PNG_INTERLACE_NONE or
   * PNG_INTERLACE_ADAM7, and the compression_type and filter_type MUST
   * currently be PNG_COMPRESSION_TYPE_BASE and PNG_FILTER_TYPE_BASE. REQUIRED
   */
  unsigned int png_color_type;
  if(channel==3) png_color_type=PNG_COLOR_TYPE_RGB;
  if(channel==4) png_color_type=PNG_COLOR_TYPE_RGB_ALPHA;

  png_set_IHDR(png_ptr, info_ptr, width, height, 8, png_color_type,
	       PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

  /* set the palette if there is one.  REQUIRED for indexed-color images */
  // palette = (png_colorp)png_malloc(png_ptr, PNG_MAX_PALETTE_LENGTH * sizeof (png_color));
  /* ... set palette colors ... */
  // png_set_PLTE(png_ptr, info_ptr, palette, PNG_MAX_PALETTE_LENGTH);
  /* You must not free palette here, because png_set_PLTE only makes a link to
     the palette that you malloced.  Wait until you are about to destroy
     the png structure. */

  /* optional significant bit chunk */
  /* if we are dealing with a grayscale image then */
  //sig_bit.gray = true_bit_depth;
  /* otherwise, if we are dealing with a color image then */
  /* png_color_8p sig_bit;

     sig_bit.red = true_red_bit_depth;
     sig_bit.green = true_green_bit_depth;
     sig_bit.blue = true_blue_bit_depth;*/
  /* if the image has an alpha channel then */
  //sig_bit.alpha = true_alpha_bit_depth;
  //   png_set_sBIT(png_ptr, info_ptr, sig_bit);

  /* Optional gamma chunk is strongly suggested if you have any guess
   * as to the correct gamma of the image.
   */
  //png_set_gAMA(png_ptr, info_ptr, gamma);

  /* Optionally write comments into the image */
  /* text_ptr[0].key = "Title";
     text_ptr[0].text = "Mona Lisa";
     text_ptr[0].compression = PNG_TEXT_COMPRESSION_NONE;
     text_ptr[1].key = "Author";
     text_ptr[1].text = "Leonardo DaVinci";
     text_ptr[1].compression = PNG_TEXT_COMPRESSION_NONE;
     text_ptr[2].key = "Description";
     text_ptr[2].text = "<long text>";
     text_ptr[2].compression = PNG_TEXT_COMPRESSION_zTXt;
     #ifdef PNG_iTXt_SUPPORTED
     text_ptr[0].lang = NULL;
     text_ptr[1].lang = NULL;
     text_ptr[2].lang = NULL;
     #endif
     png_set_text(png_ptr, info_ptr, text_ptr, 3);*/

  /* other optional chunks like cHRM, bKGD, tRNS, tIME, oFFs, pHYs, */
  /* note that if sRGB is present the gAMA and cHRM chunks must be ignored
   * on read and must be written in accordance with the sRGB profile */

  /* Write the file header information.  REQUIRED */
  png_write_info(png_ptr, info_ptr);

  /* If you want, you can write the info in two steps, in case you need to
   * write your private chunk ahead of PLTE:
   *
   *   png_write_info_before_PLTE(write_ptr, write_info_ptr);
   *   write_my_chunk();
   *   png_write_info(png_ptr, info_ptr);
   *
   * However, given the level of known- and unknown-chunk support in 1.1.0
   * and up, this should no longer be necessary.
   */

  /* Once we write out the header, the compression type on the text
   * chunks gets changed to PNG_TEXT_COMPRESSION_NONE_WR or
   * PNG_TEXT_COMPRESSION_zTXt_WR, so it doesn't get written out again
   * at the end.
   */

  /* set up the transformations you want.  Note that these are
   * all optional.  Only call them if you want them.
   */

  /* invert monochrome pixels */
  //png_set_invert_mono(png_ptr);

  /* Shift the pixels up to a legal bit depth and fill in
   * as appropriate to correctly scale the image.
   */
  //png_set_shift(png_ptr, &sig_bit);

  /* pack pixels into bytes */
  png_set_packing(png_ptr);

  /* swap location of alpha bytes from ARGB to RGBA */
  //png_set_swap_alpha(png_ptr);

  /* Get rid of filler (OR ALPHA) bytes, pack XRGB/RGBX/ARGB/RGBA into
   * RGB (4 channels -> 3 channels). The second parameter is not used.
   */
  //png_set_filler(png_ptr, 0, PNG_FILLER_BEFORE);

  /* flip BGR pixels to RGB */
  //png_set_bgr(png_ptr);

  /* swap bytes of 16-bit files to most significant byte first */
  //png_set_swap(png_ptr);

  /* swap bits of 1, 2, 4 bit packed pixel formats */
  //png_set_packswap(png_ptr);

  /* turn on interlace handling if you are not using png_write_image() */
  /* if (interlacing)
     number_passes = png_set_interlace_handling(png_ptr);
     else
     number_passes = 1;*/

  /* The easiest way to write the image (you may have a different memory
   * layout, however, so choose what fits your needs best).  You need to
   * use the first method if you aren't handling interlacing yourself.
   */
  png_uint_32 k;
  //   png_byte image[ysize][xsize*3];
  png_bytep row_pointers[height];
  for (k = 0; k < (unsigned)height; k++)
    row_pointers[k] = img + k*width*channel;

  png_write_image(png_ptr, row_pointers);

  /* You can write optional chunks like tEXt, zTXt, and tIME at the end
   * as well.  Shouldn't be necessary in 1.1.0 and up as all the public
   * chunks are supported and you can use png_set_unknown_chunks() to
   * register unknown chunks into the info structure to be written out.
   */

  /* It is REQUIRED to call this to finish writing the rest of the file */
  png_write_end(png_ptr, info_ptr);

  /* If you png_malloced a palette, free it here (don't free info_ptr->palette,
     as recommended in versions 1.0.5m and earlier of this example; if
     libpng mallocs info_ptr->palette, libpng will free it).  If you
     allocated it with malloc() instead of png_malloc(), use free() instead
     of png_free(). */


  /* Similarly, if you png_malloced any data that you passed in with
     png_set_something(), such as a hist or trans array, free it here,
     when you can be sure that libpng is through with it. */
  /*png_free(png_ptr, trans);
    trans=NULL;*/

  /* clean up after the write, and free any memory allocated */
  png_destroy_write_struct(&png_ptr, &info_ptr);

  /* close the file */
  fclose(fp);

  /* that's it */
  return (0);
}


int readPng(char *filename,
	    int *width, int *height, int *channels,
	    unsigned char **data)
{
  FILE *f;
  png_struct *png_ptr;
  png_info *info_ptr;
  png_bytep *row_pointers, row;
  int bit_depth, color_type;

  if((f=fopen(filename,"rb"))==NULL){
    printf("erreur au chargement du fichier png %s\n",filename);
    return ERR_OPENING;
  }

  /* initialiser les structures de données */

  /* initialisation de la structure png_struct */
  png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,
				   NULL, NULL, NULL);
  if(png_ptr==NULL){
    printf("erreur initialisation png_struct\n");
    fclose(f);
    return ERR_MEMORY;
  }


  /* allocation et initialisation de la structure qui */
  /* contiendra les informations sur l'image */
  info_ptr = png_create_info_struct(png_ptr);
  if(info_ptr==NULL){
    printf("erreur mémoire pour info_ptr\n");
    fclose(f);
    png_destroy_read_struct(&png_ptr, 0, 0);
    return ERR_MEMORY;
  }

  /* préparer l'adresse du code de gestion d'erreur si l'une des */
  /* fonctions de la librairie png se plante ... */ 
  if (setjmp(png_jmpbuf(png_ptr))){
    /* Free all of the memory associated with the png_ptr and info_ptr */
    png_destroy_read_struct(&png_ptr, &info_ptr, 0);
    fclose(f);
    /* If we get here, we had a problem reading the file */
    return  ERR_READING;
  }


  /* initialiser la méthode d'entrée/sortie */
  png_init_io(png_ptr, f);

  png_read_png(png_ptr, info_ptr,
	       PNG_TRANSFORM_STRIP_16 | //PNG_TRANSFORM_STRIP_ALPHA |
	       PNG_TRANSFORM_PACKING | PNG_TRANSFORM_EXPAND,
	       NULL);

  /* get some usefull information from header */
  bit_depth = png_get_bit_depth (png_ptr, info_ptr);
  color_type = png_get_color_type (png_ptr, info_ptr);

  printf("bit_depth = %d color_type = %d\n", bit_depth, color_type);

  row_pointers = png_get_rows(png_ptr, info_ptr);

  /* raz de l'image precedente */
  if(*data) delete [] *data;

  /* allocation de la nouvelle image */
  *height = png_get_image_height(png_ptr, info_ptr);
  *width = png_get_image_width(png_ptr, info_ptr);

  *channels = png_get_channels(png_ptr, info_ptr);

  cout << "nb canaux = " << *channels << endl;

  *data = new unsigned char[(*width)*(*height)*(*channels)];
  
  if(!(*data)){
    printf("erreur allocation de l'image\n");
    *width = *height = 0;
    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
    fclose(f);
    return ERR_MEMORY;    
  } 

  printf("image %dx%d\n", *width,*height);


  unsigned char *im;
  for(int y=0; y<*height; y++){
    im = *data + ((*height) -1 - y)*(*width)*(*channels) ;

    //     printf("ligne %d -", y); fflush(stdout);
    /* recupération de la ligne y */
    row = row_pointers[y];

    for(int x=0; x<(*width); x++){

      for(int c=0; c<*channels; c++){
	*im = *row;
	im++; row++;
      }

    }
  }


  png_destroy_read_struct(&png_ptr, &info_ptr, 0);

  fclose(f);

  return LOADED;

};
