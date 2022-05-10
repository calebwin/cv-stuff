#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"
#include <time.h>

// Frees an array of descriptors.
// descriptor *d: the array.
// int n: number of elements in array.
void free_descriptors(descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        free(d[i].data);
    }
    free(d);
}

// Create a feature descriptor for an index in an image.
// image im: source image.
// int i: index in image for the pixel we want to describe.
// returns: descriptor for that index.
descriptor describe_index(image im, int i)
{
    int w = 5;
    descriptor d;
    d.p.x = i%im.w;
    d.p.y = i/im.w;
    d.data = calloc(w*w*im.c, sizeof(float));
    d.n = w*w*im.c;
    int c, dx, dy;
    int count = 0;
    // If you want you can experiment with other descriptors
    // This subtracts the central value from neighbors
    // to compensate some for exposure/lighting changes.
    for(c = 0; c < im.c; ++c){
        float cval = im.data[c*im.w*im.h + i];
        for(dx = -w/2; dx < (w+1)/2; ++dx){
            for(dy = -w/2; dy < (w+1)/2; ++dy){
                float val = get_pixel(im, i%im.w+dx, i/im.w+dy, c);
                d.data[count++] = cval - val;
            }
        }
    }
    return d;
}

// Marks the spot of a point in an image.
// image im: image to mark.
// ponit p: spot to mark in the image.
void mark_spot(image im, point p)
{
    int x = p.x;
    int y = p.y;
    int i;
    for(i = -9; i < 10; ++i){
        set_pixel(im, x+i, y, 0, 1);
        set_pixel(im, x, y+i, 0, 1);
        set_pixel(im, x+i, y, 1, 0);
        set_pixel(im, x, y+i, 1, 0);
        set_pixel(im, x+i, y, 2, 1);
        set_pixel(im, x, y+i, 2, 1);
    }
}

// Marks corners denoted by an array of descriptors.
// image im: image to mark.
// descriptor *d: corners in the image.
// int n: number of descriptors to mark.
void mark_corners(image im, descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        mark_spot(im, d[i].p);
    }
}

float g(float sigma, int x) {
	return pow(M_E, -((pow(x, 2)) / (2 * pow(sigma, 2)))) / (2 * M_PI * pow(sigma, 2));
}

// Creates a 1d Gaussian filter.
// float sigma: standard deviation of Gaussian.
// returns: single row image of the filter.
image make_1d_gaussian(float sigma)
{
    int N = (int) ceil(sigma * 6.0);
    N = ((N % 2) == 0) ? N + 1 : N;
    image res = make_image(N,1,1);

    for (int x = 0; x < res.w; x++) {
        int x_relative_to_center = x - (int) floor(res.w / 2);
        float v = g(sigma, x_relative_to_center);
        set_pixel(res, x, 0, 0, v);
	}

    // Normalize
	l1_normalize(res);

    return res;
}

// Smooths an image using separable Gaussian filter.
// image im: image to smooth.
// float sigma: std dev. for Gaussian.
// returns: smoothed image.
image smooth_image(image im, float sigma)
{
    // Get filters
    image gaussian_filter_nx1 = make_1d_gaussian(sigma);
    image gaussian_filter_1xn = copy_image(gaussian_filter_nx1);
    gaussian_filter_1xn.h = gaussian_filter_1xn.w;
    gaussian_filter_1xn.w = 1;

    // Convolve
    image im2 = convolve_image(im, gaussian_filter_1xn, 1);
    image im3 = convolve_image(im2, gaussian_filter_nx1, 1);

    return im3;
}

// Calculate the structure matrix of an image.
// image im: the input image.
// float sigma: std dev. to use for weighted sum.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          third channel is IxIy.
image structure_matrix(image im, float sigma)
{
    image S = make_image(im.w, im.h, 3);
    
    image Ix_filter = make_gx_filter();
	image Iy_filter = make_gy_filter();
	image Ix = convolve_image(im, Ix_filter, 0);
	image Iy = convolve_image(im, Iy_filter, 0);
    for (int i = 0; i < S.w; i++) {
		for (int j = 0; j < S.h; j++) {
            float Ix_ij = get_pixel(Ix, i, j, 0);
            float Iy_ij = get_pixel(Iy, i, j, 0);
            set_pixel(S, i, j, 0, Ix_ij * Ix_ij);
            set_pixel(S, i, j, 1, Iy_ij * Iy_ij);
            set_pixel(S, i, j, 2, Ix_ij * Iy_ij);
        }
    }

    // image gaussian_filter = make_gaussian_filter(sigma);
    // return convolve_image(S, gaussian_filter, 1);

    return smooth_image(S, sigma);
}

// Estimate the cornerness of each pixel given a structure matrix S.
// image S: structure matrix for an image.
// returns: a response map of cornerness calculations.
image cornerness_response(image S)
{
    image R = make_image(S.w, S.h, 1);
    // We'll use formulation det(S) - alpha * trace(S)^2, alpha = .06.
    for (int i = 0; i < S.w; i++) {
		for (int j = 0; j < S.h; j++) {
            float Ixx_ij = get_pixel(S, i, j, 0);
            float Iyy_ij = get_pixel(S, i, j, 1);
            float Ixy_ij = get_pixel(S, i, j, 2);
            float Det_S_ij = (Ixx_ij * Iyy_ij) - (Ixy_ij * Ixy_ij);
            float Tr_S_ij = Ixx_ij + Iyy_ij;
            float R_ij = Det_S_ij - 0.06 * Tr_S_ij * Tr_S_ij;
            set_pixel(R, i, j, 0, R_ij);
        }
    }
    return R;
}

// Perform non-max supression on an image of feature responses.
// image im: 1-channel image of feature responses.
// int w: distance to look for larger responses.
// returns: image with only local-maxima responses within w pixels.
image nms_image(image im, int w)
{
    image r = copy_image(im);
    // for every pixel in the image:
    //     for neighbors within w:
    //         if neighbor response greater than pixel response:
    //             set response to be very low (I use -999999 [why not 0??])
    // printf("im.c=%d\n", im.c);
    for (int i = 0; i < im.w; i++) {
		for (int j = 0; j < im.h; j++) {
            int k = 0;
            for (int n_i = i - w; n_i <= i + w; n_i++) {
                for (int n_j = j - w; n_j <= j + w; n_j++) {
                    int n_k = 0;
                    // printf("get_pixel(r, i, j, k)=%f\n", get_pixel(r, i, j, k));
                    if (n_i >= 0 && n_i <= im.w - 1 && n_j >= 0 && n_j <= im.h - 1) {
                        if (get_pixel(im, n_i, n_j, n_k) > get_pixel(im, i, j, k)) {
                            // printf("i=%d, j=%d\n", i, j);
                            set_pixel(r, i, j, k, -999999);
                        }
                    }
                }
            }
        }
    }
    return r;
}

// Perform harris corner detection and extract features from the corners.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
// int *n: pointer to number of corners detected, should fill in.
// returns: array of descriptors of the corners in the image.
descriptor *harris_corner_detector(image im, float sigma, float thresh, int nms, int *n)
{
    // Calculate structure matrix
    image S = structure_matrix(im, sigma);

    // Estimate cornerness
    image R = cornerness_response(S);

    // Run NMS on the responses
    image Rnms = nms_image(R, nms);

    int count = 0;
    for (int i = 0; i < Rnms.w; i++) {
		for (int j = 0; j < Rnms.h; j++) {
            if (get_pixel(Rnms, i, j, 0) > thresh) {
                count++;
            }
        }
    }
    printf("Rnms.c=%d\n", Rnms.c);
    printf("count=%d\n", count);
    printf("thresh=%f\n", thresh);
    
    *n = count; // <- set *n equal to number of corners in image.
    descriptor *d = calloc(count, sizeof(descriptor));
    
    int di = 0;
    for (int i = 0; i < Rnms.w; i++) {
		for (int j = 0; j < Rnms.h; j++) {
            if (get_pixel(Rnms, i, j, 0) > thresh) {
                descriptor new_d = describe_index(im, (Rnms.w * j) + i);
                d[di] = new_d;
                di++;
            }
        }
    }

    free_image(S);
    free_image(R);
    free_image(Rnms);
    return d;
}

// Find and draw corners on an image.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
void detect_and_draw_corners(image im, float sigma, float thresh, int nms)
{
    int n = 0;
    descriptor *d = harris_corner_detector(im, sigma, thresh, nms, &n);
    mark_corners(im, d, n);
}
