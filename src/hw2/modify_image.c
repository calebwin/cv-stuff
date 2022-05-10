#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

/******************************** Resizing *****************************
  To resize we'll need some interpolation methods and a function to create
  a new image and fill it in with our interpolation methods.
************************************************************************/

float nn_interpolate(image im, float x, float y, int c)
{
	/***********************************************************************
	  This function performs nearest-neighbor interpolation on image "im"
	  given a floating column value "x", row value "y" and integer channel "c",
	  and returns the interpolated value.
	************************************************************************/
	return get_pixel(im, (int) round(x), (int) round(y), c);
}

image nn_resize(image im, int w, int h)
{
	/***********************************************************************
	  This function uses nearest-neighbor interpolation on image "im" to a new
	  image of size "w x h"
	************************************************************************/
	image new_im = make_image(w, h, im.c);
	float x_a = im.w * 1.0 / w;
	float x_b = -0.5 - (x_a * -0.5);
	float y_a = im.h * 1.0 / h;
	float y_b = -0.5 - (y_a * -0.5);
	for (int k = 0; k < im.c; k++)
		for (int i = 0; i < w; i++)
			for (int j = 0; j < h; j++)
				set_pixel(new_im, i, j, k, nn_interpolate(im, x_a * i + x_b, y_a * j + y_b, k));
	return new_im;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
	/***********************************************************************
		 This function performs bilinear interpolation on image "im" given
		a floating column value "x", row value "y" and integer channel "c".
		It interpolates and returns the interpolated value.
	************************************************************************/
	int x_floor = (int) floor(x);
	int y_floor = (int) floor(y);
	int x_ceil = (int) ceil(x);
	int y_ceil = (int) ceil(y);
	float V1 = get_pixel(im, x_floor, y_floor, c);
	float V2 = get_pixel(im, x_ceil, y_floor, c);
	float V3 = get_pixel(im, x_floor, y_ceil, c);
	float V4 = get_pixel(im, x_ceil, y_ceil, c);
	float d1 = x - x_floor;
	float d2 = x_ceil - x;
	float d3 = y - y_floor;
	float d4 = y_ceil - y;
	float A4 = d1 * d3;
	float A3 = d2 * d3;
	float A1 = d2 * d4;
	float A2 = d1 * d4;
	float q = V1*A1 + V2*A2 + V3*A3 + V4*A4;
	return q;
}

image bilinear_resize(image im, int w, int h)
{
	/***********************************************************************
	  This function uses bilinear interpolation on image "im" to a new image
	  of size "w x h". Algorithm is same as nearest-neighbor interpolation.
	************************************************************************/
	image new_im = make_image(w, h, im.c);
	float x_a = im.w * 1.0 / w;
	float x_b = -0.5 - (x_a * -0.5);
	float y_a = im.h * 1.0 / h;
	float y_b = -0.5 - (y_a * -0.5);
	for (int k = 0; k < im.c; k++)
		for (int i = 0; i < w; i++)
			for (int j = 0; j < h; j++)
				set_pixel(new_im, i, j, k, bilinear_interpolate(im, x_a * i + x_b, y_a * j + y_b, k));
	return new_im;
}

/********************** Filtering: Box filter ***************************
  We want to create a box filter. We will only use square box filters.
************************************************************************/

void l1_normalize(image im)
{
	/***********************************************************************
	  This function divides each value in image "im" by the sum of all the
	  values in the image and modifies the image in place.
	************************************************************************/
	float total_value = 0.0;
	for (int i = 0; i < im.w * im.h * im.c; i++) {
		total_value += im.data[i];
	}
	for (int i = 0; i < im.w * im.h * im.c; i++) {
		im.data[i] = im.data[i] / total_value;
	}
}

image make_box_filter(int w)
{
	/***********************************************************************
	  This function makes a square filter of size "w x w". Make an image of
	  width = height = w and number of channels = 1, with all entries equal
	  to 1. Then use "l1_normalize" to normalize your filter.
	************************************************************************/
	image res = make_image(w, w, 1);

	for (int i = 0; i < res.w * res.h * res.c; i++) {
		res.data[i] = 1.0;
	}

	l1_normalize(res);

	return res;
}

image convolve_image(image im, image filter, int preserve)
{
	/***********************************************************************
	  This function convolves the image "im" with the "filter". The value
	  of preserve is 1 if the number of input image channels need to be
	  preserved. Check the detailed algorithm given in the README.
	************************************************************************/
	assert(filter.c == im.c || filter.c == 1);
	
	image res = make_image(im.w, im.h, preserve ? im.c : 1);
	// int sum_across_z = im.c == filter.c &&
	int z_center_start = (preserve == 1) ? 0 : (int) floor(im.c / 2);
	int z_center_end = (preserve == 1) ? im.c : ((int) floor(im.c / 2) + 1);
	// int sum_across_z = !preserve || (filter.c == 1);
	// printf("im.c = %d, im.c/2 = %d, floor(im.c/2) = %d, preserve = %d, preserve ? 0 : 10 = %d\n", im.c, im.c/2, (int) floor(im.c/2), preserve, preserve == 1 ? 0 : 10);
	for (int z_center = z_center_start; z_center < z_center_end; z_center++) {
		for (int x_center = 0; x_center < im.w; x_center++) {
			for (int y_center = 0; y_center < im.h; y_center++) {
				float q = 0.0;
				int z_start = (preserve == 1) ? z_center : (z_center - (int) floor(im.c / 2));
				int z_end = (preserve == 1) ? (z_center + 1) : (z_center + (int) floor(im.c / 2) + 1);
				for (int z = z_start; z < z_end; z++) {
					// printf("res.h=%d, res.c=%d, z_start=%d, z_end=%d, z_center=%d, z_center_start=%d, z_center_end=%d\n", res.h, res.c, z_start, z_end, z_center, z_center_start, z_center_end);
					for (int x = x_center - (int) floor(filter.w / 2); x < (x_center + (int) floor(filter.w / 2) + 1); x++) {
						for (int y = y_center - (int) floor(filter.h / 2); y < (y_center + (int) floor(filter.h / 2) + 1); y++) {
							int z_in_filter = (filter.c == 1) ? 0 : (z - z_start);
							int x_in_filter = x - (x_center - (int) floor(filter.w / 2));
							int y_in_filter = y - (y_center - (int) floor(filter.h / 2));
							q += (
								get_pixel(im, x, y, z) *
								get_pixel(filter, x_in_filter, y_in_filter, z_in_filter)
							);
						}
					}
				}
				set_pixel(res, x_center, y_center, preserve ? z_center : 0, q);
			}
		}
	}
	return res;
}

image make_highpass_filter()
{
	/***********************************************************************
	  Create a 3x3 filter with highpass filter values using image.data[]
	************************************************************************/
	image res = make_image(3, 3, 1);
	float res_data[] = {0.0, -1.0, 0.0, -1.0, 4.0, -1.0, 0.0, -1.0, 0.0};
	memcpy(res.data, &res_data, 9 * sizeof(float));
	return res;
}

image make_sharpen_filter()
{
	/***********************************************************************
	  Create a 3x3 filter with sharpen filter values using image.data[]
	************************************************************************/
	image res = make_image(3, 3, 1);
	float res_data[] = {0.0, -1.0, 0.0, -1.0, 5.0, -1.0, 0.0, -1.0, 0.0};
	memcpy(res.data, &res_data, 9 * sizeof(float));
	return res;
}

image make_emboss_filter()
{
	/***********************************************************************
	  Create a 3x3 filter with emboss filter values using image.data[]
	************************************************************************/
	image res = make_image(3, 3, 1);
	float res_data[] = {-2.0, -1.0, 0.0, -1.0, 1.0, 1.0, 0.0, 1.0, 2.0};
	memcpy(res.data, &res_data, 9 * sizeof(float));
	return res;
}

// Question 2.3.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: Sharpen and emboss require preserving all channels while highpass doesn't because it is for edge detection where
// edges are more pronounced when the result is a single-channel image.

// Question 2.3.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: They all require clamping because the convolutions could result in values with magnitude greater than 1.

float G(float sigma, int x, int y) {
	return pow(M_E, -((pow(x, 2) + pow(y, 2)) / (2 * pow(sigma, 2)))) / (2 * M_PI * pow(sigma, 2));
}

image make_gaussian_filter(float sigma)
{
	/***********************************************************************
	  sigma: a float number for the Gaussian.
	  Create a Gaussian filter with the given sigma. Note that the kernel size
	  is the next highest odd integer from 6 x sigma. Return the Gaussian filter.
	************************************************************************/
	
	// Create filter image of the right size
	float filter_size_exact = sigma * 6.0;
	int filter_size = (int) ceil(filter_size_exact);
	filter_size = ((filter_size % 2) == 0) ? filter_size + 1 : filter_size;
	image res = make_image(filter_size, filter_size, 1);

	// Compute Gaussian distribution
	for (int x = 0; x < res.w; x++) {
		for (int y = 0; y < res.h; y++) {
			int x_relative_to_center = x - (int) floor(res.w / 2);
			int y_relative_to_center = y - (int) floor(res.h / 2);
			float v = G(sigma, x_relative_to_center, y_relative_to_center);
			set_pixel(res, x, y, 0, v);
		}
	}

	// Normalize
	l1_normalize(res);

	return res;
}

image add_image(image a, image b)
{
	/***********************************************************************
	  The input images a and image b have the same height, width, and channels.
	  Sum the given two images and return the result, which should also have
	  the same height, width, and channels as the inputs. Do necessary checks.
	************************************************************************/
	
	assert(a.w == b.w && a.h == b.h && a.c == b.c);

	image res = make_image(a.w, a.h, a.c);

	for (int i = 0; i < a.w; i++) {
		for (int j = 0; j < a.h; j++) {
			for (int k = 0; k < a.c; k++) {
				set_pixel(res, i, j, k, get_pixel(a, i, j, k) + get_pixel(b, i, j, k));
			}
		}
	}

	return res;
}

image sub_image(image a, image b)
{
	/***********************************************************************
	  The input image a and image b have the same height, width, and channels.
	  Subtract the given two images and return the result, which should have
	  the same height, width, and channels as the inputs. Do necessary checks.
	************************************************************************/
	
	assert(a.w == b.w && a.h == b.h && a.c == b.c);

	image res = make_image(a.w, a.h, a.c);

	for (int i = 0; i < a.w; i++) {
		for (int j = 0; j < a.h; j++) {
			for (int k = 0; k < a.c; k++) {
				set_pixel(res, i, j, k, get_pixel(a, i, j, k) - get_pixel(b, i, j, k));
			}
		}
	}

	return res;
}

image make_gx_filter()
{
	/***********************************************************************
	  Create a 3x3 Sobel Gx filter and return it
	************************************************************************/
	image res = make_image(3, 3, 1);
	float res_data[] = {-1.0, 0.0, 1.0, -2.0, 0.0, 2.0, -1.0, 0.0, 1.0};
	memcpy(res.data, &res_data, 9 * sizeof(float));
	return res;
}

image make_gy_filter()
{
	/***********************************************************************
	  Create a 3x3 Sobel Gy filter and return it
	************************************************************************/
	image res = make_image(3, 3, 1);
	float res_data[] = {-1.0, -2.0, -1.0, 0.0, 0.0, 0.0, 1.0, 2.0, 1.0};
	memcpy(res.data, &res_data, 9 * sizeof(float));
	return res;
}

void feature_normalize(image im)
{
	/***********************************************************************
	  Calculate minimum and maximum pixel values. Normalize the image by
	  subtracting the minimum and dividing by the max-min difference.
	************************************************************************/
	if (im.w == 0 || im.h == 0 || im.c == 0) {
		return im;
	}

	float min_value = get_pixel(im, 0, 0, 0);
	float max_value = get_pixel(im, 0, 0, 0);
	for (int i = 0; i < im.w; i++) {
		for (int j = 0; j < im.h; j++) {
			for (int k = 0; k < im.c; k++) {
				min_value = MIN(min_value, get_pixel(im, i, j, k));
				max_value = MAX(max_value, get_pixel(im, i, j, k));
			}
		}
	}
	float range = max_value - min_value;
	float multiplier = range == 0.0 ? 0.0 : 1 / range;
	for (int i = 0; i < im.w; i++) {
		for (int j = 0; j < im.h; j++) {
			for (int k = 0; k < im.c; k++) {
				float v = get_pixel(im, i, j, k);
				set_pixel(im, i, j, k, (v - min_value) * multiplier);
			}
		}
	}
}

image *sobel_image(image im)
{
	/***********************************************************************
	  Apply Sobel filter to the input image "im", get the magnitude as sobelimg[0]
	  and gradient as sobelimg[1], and return the result.
	************************************************************************/

	image Gx_filter = make_gx_filter();
	image Gy_filter = make_gy_filter();
	image Gx = convolve_image(im, Gx_filter, 0);
	image Gy = convolve_image(im, Gy_filter, 0);

	image *sobelimg = calloc(2, sizeof(image));
	sobelimg[0] = make_image(im.w, im.h, 1);
	sobelimg[1] = make_image(im.w, im.h, 1);

	for (int i = 0; i < sobelimg[0].w; i++) {
		for (int j = 0; j < sobelimg[0].h; j++) {
			float Gx_value = get_pixel(Gx, i, j, 0);
			float Gy_value = get_pixel(Gy, i, j, 0);
			set_pixel(sobelimg[0], i, j, 0, sqrt(pow(Gx_value, 2) + pow(Gy_value, 2)));
			set_pixel(sobelimg[1], i, j, 0, atan2(Gy_value, Gx_value));
		}
	}

	free(Gx_filter.data);
	free(Gy_filter.data);
	free(Gx.data);
	free(Gy.data);

	return sobelimg;
}

image colorize_sobel(image im)
{
	/***********************************************************************
	  Create a colorized version of the edges in image "im" using the
	  algorithm described in the README.
	************************************************************************/

	assert(im.c == 3);

	image *sobelimg = sobel_image(im);
	feature_normalize2(sobelimg[0]);
	feature_normalize2(sobelimg[1]);

	image res = make_image(im.w, im.h, im.c);
	for (int i = 0; i < im.w; i++) {
		for (int j = 0; j < im.h; j++) {
			float hue = get_pixel(sobelimg[1], i, j, 0);
			float saturation = get_pixel(sobelimg[0], i, j, 0);
			float value = saturation;
			set_pixel(res, i, j, 0, hue);
			set_pixel(res, i, j, 1, saturation);
			set_pixel(res, i, j, 2, value);
		}
	}

	// clamp_image(res);

	hsv_to_rgb(res);

	free(sobelimg[0].data);
	free(sobelimg[1].data);

	return res;
}

// EXTRA CREDIT: Median filter

/*
image apply_median_filter(image im, int kernel_size)
{
  return make_image(1,1,1);
}
*/

// SUPER EXTRA CREDIT: Bilateral filter

/*
image apply_bilateral_filter(image im, float sigma1, float sigma2)
{
  return make_image(1,1,1);
}
*/
