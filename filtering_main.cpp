
#include "io_bmp.h"
#include "image_comps.h"
#include<math.h>
#include<cmath>
#include<time.h>

const float PI = 3.14F;

float sinc(float sample) {
	if (sample == 0)
		return 0.4F;
	else
		return (0.4F*sinf(PI * sample * 0.4F) / (PI * sample * 0.4F));

}

float hann_win(float x, int wnd_size) {
	float sum = 0.5F * (1.0F + cosf(PI * x / float(wnd_size)));
	return 0.5F * (1.0F + cosf(PI * x / float(wnd_size)));
}


void matrix_mul(float* mat_1, float* mat_2, float* result, int filter_H) {
	float sum = 0.0F;
	for (int x = 0; x <= 2 * filter_H; x++) {
		for (int y = 0; y <= 2 * filter_H; y++) {
			result[x * (2 * filter_H + 1) + y] = mat_1[x] * mat_2[y];
			sum += result[x * (2 * filter_H + 1) + y];
			
		}

	}

	for (int i = 0; i < (2 * filter_H + 1)*(2 * filter_H + 1); i++) {

		result[i] = result[i] / sum;
	}

}

float inner_product(float* ip, int ip_stride, float* mirror_psf, int filter_extent) {
	float sum = 0.0F;
	
	for (int y = -filter_extent; y <= filter_extent; y++) {
		for (int x = -filter_extent; x <= filter_extent; x++) {
			ip[y * ip_stride + x];
			mirror_psf[y * (2 * filter_extent + 1) + x];
			sum += ip[y * ip_stride + x] * mirror_psf[y * (2 * filter_extent + 1) + x];
			
		}
	}
	return sum;
}



void my_image_comp::perform_boundary_extension()
{
  int r, c;

  // First extend upwards
  float *first_line = buf;
  for (r = 1; r <= border; r++)
	  for (c = 0; c < width; c++)
		  first_line[-r * stride + c] = first_line[(r - 1) * stride + c];//first_line[c];

  // Now extend downwards
  float *last_line = buf+(height-1)*stride;
  for (r = 1; r <= border; r++)
	  for (c = 0; c < width; c++)
		  last_line[r * stride + c] = last_line[-(r - 1) * stride + c];//[(height - 1) * -r * stride + c];    //     r * stride + c];//0;//last_line[c];

  // Now extend all rows to the left and to the right
  float *left_edge = buf-border*stride;
  float *right_edge = left_edge + width - 1;
  for (r=height+2*border; r > 0; r--, left_edge+=stride, right_edge+=stride)
    for (c=1; c <= border; c++)
      {
		left_edge[-c] =  left_edge[c];// left_edge[0];
		right_edge[c] =  right_edge[-c];// right_edge[0];
      }
}


/* ========================================================================= */
/*                              Global Functions                             */
/* ========================================================================= */

/*****************************************************************************/
/*                                apply_filter                               */
/*****************************************************************************/

void apply_filter(my_image_comp *in, my_image_comp *out, int filter_H)
{
	float gain_q1f = 0.0F;
	float gain_q1f_halfdelay = 0.0F;
	   
	 //filter length H in 1D then 2H+1 is 1D filter length

	float* q1f_n = new float[2 * filter_H + 1];   //windowd sample of streatched sinc function for even positions
	float* q1f_n_halfdelay = new float[2 * filter_H + 1];//windowd sample of streatched sinc function for odd positions

	float* rowE_colE = new float[(2 * filter_H + 1) * (2 * filter_H + 1)];
	float* rowE_colO = new float[(2 * filter_H + 1) * (2 * filter_H + 1)];
	float* rowO_colO = new float[(2 * filter_H + 1) * (2 * filter_H + 1)];
	float* rowO_colE = new float[(2 * filter_H + 1) * (2 * filter_H + 1)];

	//get the values of windowed sinc function
	for (int i = 0; i <= 2 * filter_H; i++) {

		q1f_n[i] = sinc(float(i) - float(filter_H) ) * hann_win((float(i) - float(filter_H) ), float(filter_H) + 1.0F);

	}
	
	//get the values of windowed sinc function delayed by half sample
	for (int i = 0; i <= 2 * filter_H; i++) {

		q1f_n_halfdelay[i] = sinc(float(i) - float(filter_H) - 0.5F) * hann_win((float(i) - float(filter_H) - 0.5F), float(filter_H) + 1.0F);

	}
	
	
	//multiply to get filter matrices
	matrix_mul(q1f_n, q1f_n, rowE_colE, filter_H);
	matrix_mul(q1f_n, q1f_n_halfdelay, rowE_colO, filter_H);
	matrix_mul(q1f_n_halfdelay, q1f_n_halfdelay, rowO_colO, filter_H);
	matrix_mul(q1f_n_halfdelay, q1f_n, rowO_colE, filter_H);
		

	
	//------points at origin of  the filter matrix-----//
	rowE_colE = rowE_colE + (2 * filter_H + 1) * filter_H + filter_H;  
	rowE_colO = rowE_colO + (2 * filter_H + 1) * filter_H + filter_H;

	rowO_colO = rowO_colO + (2 * filter_H + 1) * filter_H + filter_H;
	rowO_colE = rowO_colE + (2 * filter_H + 1) * filter_H + filter_H;
	//------points at origin of  the filter matrix-----//
	
	float s = 0.0F;
	

  // Perform the convolution
  for (int r = 0; r < out->height; r++) {

	  for (int c = 0; c < out->width; c++) {

		  float* op = out->buf + r * out->stride + c;
		  
		 if (r % 2 == 0 && c % 2 == 0) {
			  float* ip = in->buf + 5 * r / 2 * in->stride + 5 * c / 2;
			  *op = inner_product(ip, in->stride, rowE_colE, filter_H);

		  }
		 

		   else if (r % 2 == 0 && c % 2 != 0) {
			  float* ip = in->buf + 5 * r / 2 * in->stride + 5 * (c -1) / 2 + 2;
			  *op = inner_product(ip, in->stride, rowE_colO, filter_H);

		  }

		  
		  else if (r % 2 != 0 && c % 2 != 0) {
			  float* ip = in->buf + (5 * (r - 1) / 2 + 2) * in->stride + (5 * (c-1) / 2 + 2);
			  *op = inner_product(ip, in->stride, rowO_colO, filter_H);

		  }

		  else if (r % 2 != 0 && c % 2 == 0) {
			  float* ip = in->buf + (5 * (r-1) / 2 + 2) * in->stride + 5 * c / 2;
			  *op = inner_product(ip, in->stride, rowO_colE, filter_H);

		  }

		 
	  }

  }
   
}

/*****************************************************************************/
/*                                    main                                   */
/*****************************************************************************/

int  main(int argc, char *argv[])
{
	clock_t start_time = clock();
 
   
  if (argc != 4)
    {
      fprintf(stderr,"Usage: %s <in bmp file> <out bmp file>\n",argv[0]);
      return -1;
    }
  int filter_H = atoi(argv[3]);
  
 
  int err_code=0;
  try {
	  // Read the input image
	  bmp_in in;
	  if ((err_code = bmp_in__open(&in, argv[1])) != 0)
		  throw err_code;



	  int width = in.cols, height = in.rows;

	  int new_width = int(ceil(float(width) * 0.4F));
	  int new_height = int(ceil(float(height) * 0.4F));

	  int n, num_comps = in.num_components;
	  my_image_comp* input_comps = new my_image_comp[num_comps];

	  ///extent boundaries
	  for (n = 0; n < num_comps; n++) {
		  input_comps[n].init(height, width, filter_H); // Leave a border of filter extent size
	  }



	  int r; // Declare row index
	  io_byte* line = new io_byte[width * num_comps];
	  for (r = height - 1; r >= 0; r--) {
		  // "r" holds the true row index we are reading, since the image is
			 // stored upside down in the BMP file.
		  if ((err_code = bmp_in__get_line(&in, line)) != 0) {
			  throw err_code;
		  }

		  for (n = 0; n < num_comps; n++)
		  {
			  io_byte* src = line + n; // Points to first sample of component n

			  float* dst = input_comps[n].buf + r * input_comps[n].stride;


			  for (int c = 0; c < width; c++, src += num_comps) {

				  dst[c] = (float)* src; // The cast to type "float" is not
					   // strictly required here, since bytes can always be
					   // converted to floats without any loss of information.
			  }

		  }
	  }
  
	  bmp_in__close(&in);

      // Allocate storage for the filtered output
      my_image_comp *output_comps = new my_image_comp[num_comps];

	  for (n = 0; n < num_comps; n++) {
		  output_comps[n].init(new_height, new_width, 0); // Don't need a border for output
	  }
        

      // Process the image, all in floating point (easy)
	  for (n = 0; n < num_comps; n++) {
		  input_comps[n].perform_boundary_extension();
	  }
        
	  for (n = 0; n < num_comps; n++) {
		  apply_filter(input_comps + n, output_comps + n, filter_H);
	  }
        

      // Write the image back out again
      bmp_out out;
      if ((err_code = bmp_out__open(&out,argv[2],new_width,new_height,num_comps)) != 0)
        throw err_code;
      for (r=new_height-1; r >= 0; r--)
        { // "r" holds the true row index we are writing, since the image is
          // written upside down in BMP files.
          for (n=0; n < num_comps; n++)
            {
              io_byte *dst = line+n; // Points to first sample of component n
              float *src = output_comps[n].buf + r * output_comps[n].stride;
			  for (int c = 0; c < new_width; c++, dst += num_comps)
			  {
				  if (src[c] > 255.0F)
					  src[c] = 255;
				  else if (src[c] <0.0F)
					  src[c] = 0;

				  *dst = (io_byte)src[c];
			  } // The cast to type "io_byte" is
                      // required here, since floats cannot generally be
                      // converted to bytes without loss of information.  The
                      // compiler will warn you of this if you remove the cast.
                      // There is in fact not the best way to do the
                      // conversion.  You should fix it up in the lab.
            }
          bmp_out__put_line(&out,line);
        }
      bmp_out__close(&out);
      delete[] line;
      delete[] input_comps;
      delete[] output_comps;
    }
  catch (int exc) {
      if (exc == IO_ERR_NO_FILE)
        fprintf(stderr,"Cannot open supplied input or output file.\n");
      else if (exc == IO_ERR_FILE_HEADER)
        fprintf(stderr,"Error encountered while parsing BMP file header.\n");
      else if (exc == IO_ERR_UNSUPPORTED)
        fprintf(stderr,"Input uses an unsupported BMP file format.\n  Current "
                "simple example supports only 8-bit and 24-bit data.\n");
      else if (exc == IO_ERR_FILE_TRUNC)
        fprintf(stderr,"Input or output file truncated unexpectedly.\n");
      else if (exc == IO_ERR_FILE_NOT_OPEN)
        fprintf(stderr,"Trying to access a file which is not open!(?)\n");
      return -1;
    }
  clock_t end_time = clock();
  float seconds = (float((end_time - start_time))) / CLOCKS_PER_SEC;
  printf("Computation Time = %f\n", seconds);
  return 0;
}
