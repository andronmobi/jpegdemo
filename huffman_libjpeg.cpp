/*
 * @file huffman_libjpeg.cpp
 *
 * @author
 * @date Mars 0x2015
 */

#include <stdio.h>
#include <stdint.h>

#include <jpeglib.h>
#include <jpegint.h>
#include <jerror.h>
#include <setjmp.h>

using namespace std;

static const char* filename = "testorig.jpg";

/*
 * Expanded entropy decoder object for Huffman decoding.
 *
 * The savable_state subrecord contains fields that change within an MCU,
 * but must not be updated permanently until we complete the MCU.
 */

struct my_error_mgr {
    struct jpeg_error_mgr pub;	/* "public" fields */

    jmp_buf setjmp_buffer;	/* for return to caller */
};

typedef struct my_error_mgr * my_error_ptr;

/*
 * Here's the routine that will replace the standard error_exit method:
 */

    METHODDEF(void)
my_error_exit (j_common_ptr cinfo)
{
    /* cinfo->err really points to a my_error_mgr struct, so coerce pointer */
    my_error_ptr myerr = (my_error_ptr) cinfo->err;

    /* Always display the message. */
    /* We could postpone this until after returning, if we chose. */
    (*cinfo->err->output_message) (cinfo);

    /* Return control to the setjmp point */
    longjmp(myerr->setjmp_buffer, 1);
}

void make_d_derived_tbl(jpeg_decompress_struct* cinfo, JHUFF_TBL* htbl) {
    int p, i, l, si, numsymbols;
    int lookbits, ctr;
    char huffsize[257];
    unsigned int huffcode[257];
    unsigned int code;

    /* Figure C.1: make table of Huffman code length for each symbol */
    p = 0;
    for (l = 1; l <= 16; l++) {
        i = (int) htbl->bits[l];
        if (i < 0 || p + i > 256)	/* protect against table overrun */
            printf("error JERR_BAD_HUFF_TABLE\n");
        while (i--) {
            printf("huffsize[%d]=%d ", p, l);
            huffsize[p++] = (char) l;
        }
    }
    huffsize[p] = 0;
    numsymbols = p;
    printf("\n");

    /* Figure C.2: generate the codes themselves */
    /* We also validate that the counts represent a legal Huffman code tree. */

    code = 0;
    si = huffsize[0];
    p = 0;
    while (huffsize[p]) {
        while (((int) huffsize[p]) == si) {
            printf("huffcode[%d]=%d ", p, code);
            huffcode[p++] = code;
            code++;
        }
        /* code is now 1 more than the last code used for codelength si; but
         * it must still fit in si bits, since no code is allowed to be all ones.
         */
        if (((INT32) code) >= (((INT32) 1) << si))
            printf("error JERR_BAD_HUFF_TABLE\n");
        code <<= 1;
        si++;
    }
    printf("\n");
}

// ***********************************************************

int main (int argc, char *argv[]) {

    printf("Huffman tables by libjpeg\n");

    /* This struct contains the JPEG decompression parameters and pointers to
     * working space (which is allocated as needed by the JPEG library).
     */
    struct jpeg_decompress_struct cinfo;
    /* We use our private extension JPEG error handler.
     * Note that this struct must live as long as the main JPEG parameter
     * struct, to avoid dangling-pointer problems.
     */
    struct my_error_mgr jerr;
    /* More stuff */
    FILE * infile;		/* source file */
    JSAMPARRAY buffer;		/* Output row buffer */
    int row_stride;		/* physical row width in output buffer */

    /* In this example we want to open the input file before doing anything else,
     * so that the setjmp() error recovery below can assume the file is open.
     * VERY IMPORTANT: use "b" option to fopen() if you are on a machine that
     * requires it in order to read binary files.
     */

    if ((infile = fopen(filename, "rb")) == NULL) {
        fprintf(stderr, "can't open %s\n", filename);
        return 0;
    }

    /* Step 1: allocate and initialize JPEG decompression object */

    /* We set up the normal JPEG error routines, then override error_exit. */
    cinfo.err = jpeg_std_error(&jerr.pub);
    jerr.pub.error_exit = my_error_exit;
    /* Establish the setjmp return context for my_error_exit to use. */
    if (setjmp(jerr.setjmp_buffer)) {
        /* If we get here, the JPEG code has signaled an error.
         * We need to clean up the JPEG object, close the input file, and return.
         */
        jpeg_destroy_decompress(&cinfo);
        fclose(infile);
        return 0;
    }
    /* Now we can initialize the JPEG decompression object. */
    jpeg_create_decompress(&cinfo);

    /* Step 2: specify data source (eg, a file) */

    jpeg_stdio_src(&cinfo, infile);

    /* Step 3: read file parameters with jpeg_read_header() */

    (void) jpeg_read_header(&cinfo, TRUE);
    /* We can ignore the return value from jpeg_read_header since
     *   (a) suspension is not possible with the stdio data source, and
     *   (b) we passed TRUE to reject a tables-only JPEG file as an error.
     * See libjpeg.doc for more info.
     */

    for (int i = 0; i < NUM_HUFF_TBLS; i++) {
        printf("dc_huff_table[%i] = %p\n", i, cinfo.dc_huff_tbl_ptrs[i]);
        if (cinfo.dc_huff_tbl_ptrs[i] != NULL) {
            for (int j = 1; j < 17; j++) {
                printf("%d ", cinfo.dc_huff_tbl_ptrs[i]->bits[j]);
            }
            printf("\n");
            for (int j = 0; j < 256; j++) {
                printf("%d ", cinfo.dc_huff_tbl_ptrs[i]->huffval[j]);
            }
            printf("\n");
            make_d_derived_tbl(&cinfo, cinfo.dc_huff_tbl_ptrs[i]);
        }
        printf("ac_huff_table[%i] = %p\n", i, cinfo.ac_huff_tbl_ptrs[i]);
        if (cinfo.ac_huff_tbl_ptrs[i] != NULL) {
            for (int j = 1; j < 17; j++) {
                printf("%d ", cinfo.ac_huff_tbl_ptrs[i]->bits[j]);
            }
            printf("\n");
            for (int j = 0; j < 256; j++) {
                printf("%d ", cinfo.ac_huff_tbl_ptrs[i]->huffval[j]);
            }
            printf("\n");
            make_d_derived_tbl(&cinfo, cinfo.ac_huff_tbl_ptrs[i]);
        }
    }

    /* Step 4: set parameters for decompression */

    /* In this example, we don't need to change any of the defaults set by
     * jpeg_read_header(), so we do nothing here.
     */

    /* Step 5: Start decompressor */

    (void) jpeg_start_decompress(&cinfo);
    /* We can ignore the return value since suspension is not possible
     * with the stdio data source.
     */

    /* We may need to do some setup of our own at this point before reading
     * the data.  After jpeg_start_decompress() we have the correct scaled
     * output image dimensions available, as well as the output colormap
     * if we asked for color quantization.
     * In this example, we need to make an output work buffer of the right size.
     */ 
    /* JSAMPLEs per row in output buffer */
    row_stride = cinfo.output_width * cinfo.output_components;
    /* Make a one-row-high sample array that will go away when done with image */
    buffer = (*cinfo.mem->alloc_sarray)
        ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);

    /* Step 6: while (scan lines remain to be read) */
    /*           jpeg_read_scanlines(...); */

    /* Here we use the library's state variable cinfo.output_scanline as the
     * loop counter, so that we don't have to keep track ourselves.
     */
    while (cinfo.output_scanline < cinfo.output_height) {
        /* jpeg_read_scanlines expects an array of pointers to scanlines.
         * Here the array is only one element long, but you could ask for
         * more than one scanline at a time if that's more convenient.
         */
        (void) jpeg_read_scanlines(&cinfo, buffer, 1);
        /* Assume put_scanline_someplace wants a pointer and sample count. */
        //put_scanline_someplace(buffer[0], row_stride); TODO
    }

    /* Step 7: Finish decompression */

    (void) jpeg_finish_decompress(&cinfo);
    /* We can ignore the return value since suspension is not possible
     * with the stdio data source.
     */

    /* Step 8: Release JPEG decompression object */

    /* This is an important step since it will release a good deal of memory. */
    jpeg_destroy_decompress(&cinfo);

    /* After finish_decompress, we can close the input file.
     * Here we postpone it until after no more JPEG errors are possible,
     * so as to simplify the setjmp error logic above.  (Actually, I don't
     * think that jpeg_destroy can do an error exit, but why assume anything...)
     */
    fclose(infile);

    return 0;
}

