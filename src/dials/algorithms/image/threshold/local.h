/*
 * local.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_IMAGE_THRESHOLD_LOCAL_H
#define DIALS_ALGORITHMS_IMAGE_THRESHOLD_LOCAL_H

#include <cmath>
#include <vector>
#include <iostream>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/ref_reductions.h>
#include <dials/error.h>
#include <dials/algorithms/image/filter/mean_and_variance.h>
#include <dials/algorithms/image/filter/index_of_dispersion_filter.h>
#include <dials/algorithms/image/filter/distance.h>

namespace dials { namespace algorithms {

  /**
   * Threshold the image using the niblack method.
   *
   * pixel > mean + n_sigma * sdev ? object : background
   *
   * @param image The image to threshold
   * @param size The size of the local area
   * @param n_sigma The number of standard deviations
   * @returns The thresholded image
   */
  template <typename FloatType>
  af::versa<bool, af::c_grid<2> > niblack(
    const af::const_ref<FloatType, af::c_grid<2> > &image,
    int2 size,
    double n_sigma) {
    // Check the input
    DIALS_ASSERT(n_sigma >= 0);

    // Calculate the mean and variance filtered images
    MeanAndVarianceFilter<FloatType> filter(image, size);
    af::versa<FloatType, af::c_grid<2> > mean = filter.mean();
    af::versa<FloatType, af::c_grid<2> > var = filter.sample_variance();

    // Assign the pixels to object and background
    af::versa<bool, af::c_grid<2> > result(image.accessor(),
                                           af::init_functor_null<bool>());
    for (std::size_t i = 0; i < var.size(); ++i) {
      result[i] = image[i] > mean[i] + n_sigma * std::sqrt(var[i]) ? 1 : 0;
    }

    // Return the thresholded image
    return result;
  }

  /**
   * Threshold the image using the sauvola method.
   *
   * pixel > mean * (1 + k * (sdev / (r - 1)))
   *
   * @param image The image to threshold
   * @param size The size of the local area
   * @param k A parameter
   * @param r A parameter
   * @returns The thresholded image
   */
  template <typename FloatType>
  af::versa<bool, af::c_grid<2> > sauvola(
    const af::const_ref<FloatType, af::c_grid<2> > &image,
    int2 size,
    double k,
    double r) {
    // Check the input
    DIALS_ASSERT(k >= 0 && r > 1);

    // Calculate the mean and variance filtered images
    MeanAndVarianceFilter<FloatType> filter(image, size);
    af::versa<FloatType, af::c_grid<2> > mean = filter.mean();
    af::versa<FloatType, af::c_grid<2> > var = filter.sample_variance();

    // Assign the pixels to object and background
    af::versa<bool, af::c_grid<2> > result(image.accessor(),
                                           af::init_functor_null<bool>());
    for (std::size_t i = 0; i < var.size(); ++i) {
      result[i] = image[i] > mean[i] * (1.0 + k * (std::sqrt(var[i]) / r - 1)) ? 1 : 0;
    }

    // Return the thresholded image
    return result;
  }

  /**
   * Threshold the image using a index_of_dispersion filter. Essentially a test for
   * objects within a poisson distribution.
   *
   * var/mean > 1.0 + n_sigma * sqrt(2 / (n - 1)) ? object : background
   *
   * @param image The image to threshold
   * @param size The size of the local window
   * @param n_sigma The number of standard deviations.
   */
  template <typename FloatType>
  af::versa<bool, af::c_grid<2> > index_of_dispersion(
    const af::const_ref<FloatType, af::c_grid<2> > &image,
    int2 size,
    double n_sigma) {
    // Check the input
    DIALS_ASSERT(n_sigma >= 0);

    // Calculate the index_of_dispersion filtered image
    IndexOfDispersionFilter<FloatType> filter(image, size);
    af::versa<FloatType, af::c_grid<2> > index_of_dispersion_image =
      filter.index_of_dispersion();

    // Calculate the bound
    std::size_t n = (2 * size[0] + 1) * (2 * size[1] + 1);
    DIALS_ASSERT(n > 1);
    FloatType bound = 1.0 + n_sigma * std::sqrt(2.0 / (n - 1));

    // Assign pixels to object or background
    af::versa<bool, af::c_grid<2> > result(image.accessor(),
                                           af::init_functor_null<bool>());
    for (std::size_t i = 0; i < image.size(); ++i) {
      result[i] = (index_of_dispersion_image[i] > bound) ? 1 : 0;
    }

    // Return thresholded image
    return result;
  }

  /**
   * Threshold the image using a index_of_dispersion filter. Essentially a test for
   * objects within a poisson distribution.
   *
   * var/mean > 1.0 + n_sigma * sqrt(2 / (n - 1)) ? object : background
   *
   * @param image The image to threshold
   * @param mask The mask to use
   * @param size The size of the local window
   * @param min_count The minimum counts for a point to be valid
   * @param n_sigma The number of standard deviations.
   */
  template <typename FloatType>
  af::versa<bool, af::c_grid<2> > index_of_dispersion_masked(
    const af::const_ref<FloatType, af::c_grid<2> > &image,
    const af::const_ref<bool, af::c_grid<2> > &mask,
    int2 size,
    int min_count,
    double n_sigma) {
    // Check the input
    DIALS_ASSERT(n_sigma >= 0);
    DIALS_ASSERT(min_count > 1);

    // Copy the mask into a temp variable
    af::versa<int, af::c_grid<2> > temp(mask.accessor());
    for (std::size_t i = 0; i < temp.size(); ++i) {
      temp[i] = mask[i] ? 1 : 0;
    }

    // Calculate the masked index_of_dispersion filtered image
    IndexOfDispersionFilterMasked<FloatType> filter(
      image, temp.const_ref(), size, min_count);
    af::versa<FloatType, af::c_grid<2> > index_of_dispersion_image =
      filter.index_of_dispersion();
    af::versa<int, af::c_grid<2> > count = filter.count();
    temp = filter.mask();

    // Assign pixels to object or background
    af::versa<bool, af::c_grid<2> > result(image.accessor(), false);
    for (std::size_t i = 0; i < image.size(); ++i) {
      if (temp[i]) {
        FloatType bound = 1.0 + n_sigma * std::sqrt(2.0 / (count[i] - 1));
        result[i] = (index_of_dispersion_image[i] > bound) ? 1 : 0;
      }
    }

    // Return thresholded image
    return result;
  }

  /**
   * Threshold the image using a gain filter. Same as the index_of_dispersion filter but
   * using a gain map for the calculation
   *
   * var/mean > g + n_sigma * g * sqrt(2 / (n - 1)) ? object : background
   *
   * @param image The image to threshold
   * @param mask The mask to use
   * @param gain The gain map
   * @param size The size of the local window
   * @param min_count The minimum counts for a point to be valid
   * @param n_sigma The number of standard deviations.
   */
  template <typename FloatType>
  af::versa<bool, af::c_grid<2> > gain(
    const af::const_ref<FloatType, af::c_grid<2> > &image,
    const af::const_ref<bool, af::c_grid<2> > &mask,
    const af::const_ref<FloatType, af::c_grid<2> > &gain,
    int2 size,
    int min_count,
    double n_sigma) {
    // Check the input
    DIALS_ASSERT(n_sigma >= 0);
    DIALS_ASSERT(min_count > 1);

    // Copy the mask into a temp variable
    af::versa<int, af::c_grid<2> > temp(mask.accessor());
    for (std::size_t i = 0; i < temp.size(); ++i) {
      temp[i] = mask[i] ? 1 : 0;
    }

    // Calculate the masked index_of_dispersion filtered image
    IndexOfDispersionFilterMasked<FloatType> filter(
      image, temp.const_ref(), size, min_count);
    af::versa<FloatType, af::c_grid<2> > index_of_dispersion_image =
      filter.index_of_dispersion();
    af::versa<int, af::c_grid<2> > count = filter.count();
    temp = filter.mask();

    // Assign pixels to object or background
    af::versa<bool, af::c_grid<2> > result(image.accessor(), false);
    for (std::size_t i = 0; i < image.size(); ++i) {
      if (temp[i]) {
        FloatType bound = gain[i] + n_sigma * gain[i] * std::sqrt(2.0 / (count[i] - 1));
        result[i] = (index_of_dispersion_image[i] > bound) ? 1 : 0;
      }
    }

    // Return thresholded image
    return result;
  }

  /**
   * Threshold the image as in xds. Same as the index_of_dispersion filter but
   * using a gain map for the calculation
   *
   * var/mean > g + n_sigma * g * sqrt(2 / (n - 1)) &&
   * pixel > mean + sqrt(mean) ? object : background
   *
   * @param image The image to threshold
   * @param mask The mask to use
   * @param size The size of the local window
   * @param nsig_b The background threshold.
   * @param nsig_s The strong pixel threshold
   * @param min_count The minimum number of pixels in the local area
   */
  template <typename FloatType>
  af::versa<bool, af::c_grid<2> > dispersion(
    const af::const_ref<FloatType, af::c_grid<2> > &image,
    const af::const_ref<bool, af::c_grid<2> > &mask,
    int2 size,
    double nsig_b,
    double nsig_s,
    int min_count) {
    // Check the input
    DIALS_ASSERT(nsig_b >= 0 && nsig_s >= 0);

    // Copy the mask into a temp variable
    af::versa<int, af::c_grid<2> > temp(mask.accessor());
    for (std::size_t i = 0; i < temp.size(); ++i) {
      temp[i] = mask[i] ? 1 : 0;
    }

    // Calculate the masked index_of_dispersion filtered image
    IndexOfDispersionFilterMasked<FloatType> filter(
      image, temp.const_ref(), size, min_count);
    af::versa<FloatType, af::c_grid<2> > index_of_dispersion_image =
      filter.index_of_dispersion();
    af::versa<FloatType, af::c_grid<2> > mean = filter.mean();
    af::versa<int, af::c_grid<2> > count = filter.count();
    temp = filter.mask();

    // Assign pixels to object or background
    af::versa<bool, af::c_grid<2> > result(image.accessor(), false);
    for (std::size_t i = 0; i < image.size(); ++i) {
      if (temp[i]) {
        FloatType bnd_b = 1.0 + nsig_b * std::sqrt(2.0 / (count[i] - 1));
        FloatType bnd_s = mean[i] + nsig_s * std::sqrt(mean[i]);
        result[i] = (index_of_dispersion_image[i] > bnd_b && image[i] > bnd_s) ? 1 : 0;
      }
    }

    // Return thresholded image
    return result;
  }

  /**
   * Threshold the image as in xds. Same as the index_of_dispersion filter but
   * using a gain map for the calculation
   *
   * var/mean > g + n_sigma * g * sqrt(2 / (n - 1)) &&
   * pixel > mean + sqrt(gain * mean) ? object : background
   *
   * @param image The image to threshold
   * @param mask The mask to use
   * @param gain The gain map
   * @param size The size of the local window
   * @param nsig_b The background threshold.
   * @param nsig_s The strong pixel threshold
   * @param min_count The minimum number of pixels in the local area
   */
  template <typename FloatType>
  af::versa<bool, af::c_grid<2> > dispersion_w_gain(
    const af::const_ref<FloatType, af::c_grid<2> > &image,
    const af::const_ref<bool, af::c_grid<2> > &mask,
    const af::const_ref<FloatType, af::c_grid<2> > &gain,
    int2 size,
    double nsig_b,
    double nsig_s,
    int min_count) {
    // Check the input
    DIALS_ASSERT(nsig_b >= 0 && nsig_s >= 0);

    // Copy the mask into a temp variable
    af::versa<int, af::c_grid<2> > temp(mask.accessor());
    for (std::size_t i = 0; i < temp.size(); ++i) {
      temp[i] = mask[i] ? 1 : 0;
    }

    // Calculate the masked index_of_dispersion filtered image
    IndexOfDispersionFilterMasked<FloatType> filter(
      image, temp.const_ref(), size, min_count);
    af::versa<FloatType, af::c_grid<2> > index_of_dispersion_image =
      filter.index_of_dispersion();
    af::versa<FloatType, af::c_grid<2> > mean = filter.mean();
    af::versa<int, af::c_grid<2> > count = filter.count();
    temp = filter.mask();

    // Assign pixels to object or background
    af::versa<bool, af::c_grid<2> > result(image.accessor(), false);
    for (std::size_t i = 0; i < image.size(); ++i) {
      if (temp[i]) {
        FloatType bnd_b = gain[i] + nsig_b * gain[i] * std::sqrt(2.0 / (count[i] - 1));
        FloatType bnd_s = mean[i] + nsig_s * std::sqrt(gain[i] * mean[i]);
        result[i] = (index_of_dispersion_image[i] > bnd_b && image[i] > bnd_s) ? 1 : 0;
      }
    }

    // Return thresholded image
    return result;
  }

  /**
   * A class to compute the threshold using index of dispersion
   */
  class DispersionThreshold {
  public:
    /**
     * Enable more efficient memory usage by putting components required for the
     * summed area table closer together in memory
     */
    template <typename T>
    struct Data {
      int m;
      T x;
      T y;
    };

    DispersionThreshold(int2 image_size,
                        int2 kernel_size,
                        double nsig_b,
                        double nsig_s,
                        double threshold,
                        int min_count)
        : image_size_(image_size),
          kernel_size_(kernel_size),
          nsig_b_(nsig_b),
          nsig_s_(nsig_s),
          threshold_(threshold),
          min_count_(min_count) {
      // Check the input
      DIALS_ASSERT(threshold_ >= 0);
      DIALS_ASSERT(nsig_b >= 0 && nsig_s >= 0);
      DIALS_ASSERT(image_size.all_gt(0));
      DIALS_ASSERT(kernel_size.all_gt(0));

      // Ensure the min counts are valid
      std::size_t num_kernel = (2 * kernel_size[0] + 1) * (2 * kernel_size[1] + 1);
      if (min_count_ <= 0) {
        min_count_ = num_kernel;
      } else {
        DIALS_ASSERT(min_count_ <= num_kernel && min_count_ > 1);
      }

      // Allocate the buffer
      std::size_t element_size = sizeof(Data<double>);
      buffer_.resize(element_size * image_size[0] * image_size[1]);
    }

    /**
     * Compute the summed area tables for the mask, src and src^2.
     * @param src The input array
     * @param mask The mask array
     */
    template <typename T>
    void compute_sat(af::ref<Data<T> > table,
                     const af::const_ref<T, af::c_grid<2> > &src,
                     const af::const_ref<bool, af::c_grid<2> > &mask) {
      // Largest value to consider
      const T BIG = (1 << 24);  // About 16m counts

      // Get the size of the image
      std::size_t ysize = src.accessor()[0];
      std::size_t xsize = src.accessor()[1];

      // Create the summed area table
      for (std::size_t j = 0, k = 0; j < ysize; ++j) {
        int m = 0;
        T x = 0;
        T y = 0;
        for (std::size_t i = 0; i < xsize; ++i, ++k) {
          int mm = (mask[k] && src[k] < BIG) ? 1 : 0;
          m += mm;
          x += mm * src[k];
          y += mm * src[k] * src[k];
          if (j == 0) {
            table[k].m = m;
            table[k].x = x;
            table[k].y = y;
          } else {
            table[k].m = table[k - xsize].m + m;
            table[k].x = table[k - xsize].x + x;
            table[k].y = table[k - xsize].y + y;
          }
        }
      }
    }

    /**
     * Compute the threshold
     * @param src - The input array
     * @param mask - The mask array
     * @param dst The output array
     */
    template <typename T>
    void compute_threshold(af::ref<Data<T> > table,
                           const af::const_ref<T, af::c_grid<2> > &src,
                           const af::const_ref<bool, af::c_grid<2> > &mask,
                           af::ref<bool, af::c_grid<2> > dst) {
      // Get the size of the image
      std::size_t ysize = src.accessor()[0];
      std::size_t xsize = src.accessor()[1];

      // The kernel size
      int kxsize = kernel_size_[1];
      int kysize = kernel_size_[0];

      // Calculate the local mean at every point
      for (std::size_t j = 0, k = 0; j < ysize; ++j) {
        for (std::size_t i = 0; i < xsize; ++i, ++k) {
          int i0 = i - kxsize - 1, i1 = i + kxsize;
          int j0 = j - kysize - 1, j1 = j + kysize;
          i1 = i1 < xsize ? i1 : xsize - 1;
          j1 = j1 < ysize ? j1 : ysize - 1;
          int k0 = j0 * xsize;
          int k1 = j1 * xsize;

          // Compute the number of points valid in the local area,
          // the sum of the pixel values and the sum of the squared pixel
          // values.
          double m = 0;
          double x = 0;
          double y = 0;
          if (i0 >= 0 && j0 >= 0) {
            const Data<T> &d00 = table[k0 + i0];
            const Data<T> &d10 = table[k1 + i0];
            const Data<T> &d01 = table[k0 + i1];
            m += d00.m - (d10.m + d01.m);
            x += d00.x - (d10.x + d01.x);
            y += d00.y - (d10.y + d01.y);
          } else if (i0 >= 0) {
            const Data<T> &d10 = table[k1 + i0];
            m -= d10.m;
            x -= d10.x;
            y -= d10.y;
          } else if (j0 >= 0) {
            const Data<T> &d01 = table[k0 + i1];
            m -= d01.m;
            x -= d01.x;
            y -= d01.y;
          }
          const Data<T> &d11 = table[k1 + i1];
          m += d11.m;
          x += d11.x;
          y += d11.y;

          // Compute the thresholds
          dst[k] = false;
          if (mask[k] && m >= min_count_ && x >= 0 && src[k] > threshold_) {
            double a = m * y - x * x - x * (m - 1);
            double b = m * src[k] - x;
            double c = x * nsig_b_ * std::sqrt(2 * (m - 1));
            double d = nsig_s_ * std::sqrt(x * m);
            dst[k] = a > c && b > d;
          }
        }
      }
    }

    /**
     * Compute the threshold
     * @param src - The input array
     * @param mask - The mask array
     * @param gain - The gain array
     * @param dst The output array
     */
    template <typename T>
    void compute_threshold(af::ref<Data<T> > table,
                           const af::const_ref<T, af::c_grid<2> > &src,
                           const af::const_ref<bool, af::c_grid<2> > &mask,
                           const af::const_ref<double, af::c_grid<2> > &gain,
                           af::ref<bool, af::c_grid<2> > dst) {
      // Get the size of the image
      std::size_t ysize = src.accessor()[0];
      std::size_t xsize = src.accessor()[1];

      // The kernel size
      int kxsize = kernel_size_[1];
      int kysize = kernel_size_[0];

      // Calculate the local mean at every point
      for (std::size_t j = 0, k = 0; j < ysize; ++j) {
        for (std::size_t i = 0; i < xsize; ++i, ++k) {
          int i0 = i - kxsize - 1, i1 = i + kxsize;
          int j0 = j - kysize - 1, j1 = j + kysize;
          i1 = i1 < xsize ? i1 : xsize - 1;
          j1 = j1 < ysize ? j1 : ysize - 1;
          int k0 = j0 * xsize;
          int k1 = j1 * xsize;

          // Compute the number of points valid in the local area,
          // the sum of the pixel values and the num of the squared pixel
          // values.
          double m = 0;
          double x = 0;
          double y = 0;
          if (i0 >= 0 && j0 >= 0) {
            const Data<T> &d00 = table[k0 + i0];
            const Data<T> &d10 = table[k1 + i0];
            const Data<T> &d01 = table[k0 + i1];
            m += d00.m - (d10.m + d01.m);
            x += d00.x - (d10.x + d01.x);
            y += d00.y - (d10.y + d01.y);
          } else if (i0 >= 0) {
            const Data<T> &d10 = table[k1 + i0];
            m -= d10.m;
            x -= d10.x;
            y -= d10.y;
          } else if (j0 >= 0) {
            const Data<T> &d01 = table[k0 + i1];
            m -= d01.m;
            x -= d01.x;
            y -= d01.y;
          }
          const Data<T> &d11 = table[k1 + i1];
          m += d11.m;
          x += d11.x;
          y += d11.y;

          // Compute the thresholds
          dst[k] = false;
          if (mask[k] && m >= min_count_ && x >= 0 && src[k] > threshold_) {
            double a = m * y - x * x;
            double b = m * src[k] - x;
            double c = gain[k] * x * (m - 1 + nsig_b_ * std::sqrt(2 * (m - 1)));
            double d = nsig_s_ * std::sqrt(gain[k] * x * m);
            dst[k] = a > c && b > d;
          }
        }
      }
    }

    /**
     * Compute the threshold for the given image and mask.
     * @param src - The input image array.
     * @param mask - The mask array.
     * @param dst - The destination array.
     */
    template <typename T>
    void threshold(const af::const_ref<T, af::c_grid<2> > &src,
                   const af::const_ref<bool, af::c_grid<2> > &mask,
                   af::ref<bool, af::c_grid<2> > dst) {
      // check the input
      DIALS_ASSERT(src.accessor().all_eq(image_size_));
      DIALS_ASSERT(src.accessor().all_eq(mask.accessor()));
      DIALS_ASSERT(src.accessor().all_eq(dst.accessor()));

      // Get the table
      DIALS_ASSERT(sizeof(T) <= sizeof(double));

      // Cast the buffer to the table type
      af::ref<Data<T> > table(reinterpret_cast<Data<T> *>(&buffer_[0]), buffer_.size());

      // compute the summed area table
      compute_sat(table, src, mask);

      // Compute the image threshold
      compute_threshold(table, src, mask, dst);
    }

    /**
     * Compute the threshold for the given image and mask.
     * @param src - The input image array.
     * @param mask - The mask array.
     * @param gain - The gain array
     * @param dst - The destination array.
     */
    template <typename T>
    void threshold_w_gain(const af::const_ref<T, af::c_grid<2> > &src,
                          const af::const_ref<bool, af::c_grid<2> > &mask,
                          const af::const_ref<double, af::c_grid<2> > &gain,
                          af::ref<bool, af::c_grid<2> > dst) {
      // check the input
      DIALS_ASSERT(src.accessor().all_eq(image_size_));
      DIALS_ASSERT(src.accessor().all_eq(mask.accessor()));
      DIALS_ASSERT(src.accessor().all_eq(gain.accessor()));
      DIALS_ASSERT(src.accessor().all_eq(dst.accessor()));

      // Get the table
      DIALS_ASSERT(sizeof(T) <= sizeof(double));

      // Cast the buffer to the table type
      af::ref<Data<T> > table((Data<T> *)&buffer_[0], buffer_.size());

      // compute the summed area table
      compute_sat(table, src, mask);

      // Compute the image threshold
      compute_threshold(table, src, mask, gain, dst);
    }

  private:
    int2 image_size_;
    int2 kernel_size_;
    double nsig_b_;
    double nsig_s_;
    double threshold_;
    int min_count_;
    std::vector<char> buffer_;
  };

  /**
   * A class to help debug spot finding by exposing the results of various bits
   * of processing.
   */
  class DispersionThresholdDebug {
  public:
    /**
     * Do the processing.
     * @param image The image array
     * @param mask The mask array
     * @param size The size of the local window
     * @param nsig_b The background threshold.
     * @param nsig_s The strong pixel threshold
     * @param threshold The global threshold value
     * @param min_count The minimum number of pixels in the local area
     */
    DispersionThresholdDebug(const af::const_ref<double, af::c_grid<2> > &image,
                             const af::const_ref<bool, af::c_grid<2> > &mask,
                             int2 size,
                             double nsig_b,
                             double nsig_s,
                             double threshold,
                             int min_count) {
      af::versa<double, af::c_grid<2> > gain(image.accessor(), 1.0);
      init(image, mask, gain.const_ref(), size, nsig_b, nsig_s, threshold, min_count);
    }

    /**
     * Do the processing.
     * @param image The image array
     * @param mask The mask array
     * @param size The size of the local window
     * @param nsig_b The background threshold.
     * @param nsig_s The strong pixel threshold
     * @param threshold The global threshold value
     * @param min_count The minimum number of pixels in the local area
     */
    DispersionThresholdDebug(const af::const_ref<double, af::c_grid<2> > &image,
                             const af::const_ref<bool, af::c_grid<2> > &mask,
                             const af::const_ref<double, af::c_grid<2> > &gain,
                             int2 size,
                             double nsig_b,
                             double nsig_s,
                             double threshold,
                             int min_count) {
      init(image, mask, gain, size, nsig_b, nsig_s, threshold, min_count);
    }

    /** @returns The mean map */
    af::versa<double, af::c_grid<2> > mean() const {
      return mean_;
    }

    /** @returns The variance map. */
    af::versa<double, af::c_grid<2> > variance() const {
      return variance_;
    }

    /** @returns The index of dispersion map */
    af::versa<double, af::c_grid<2> > index_of_dispersion() const {
      return cv_;
    }

    /** @returns The thresholded index of dispersion mask */
    af::versa<bool, af::c_grid<2> > cv_mask() const {
      return cv_mask_;
    }

    /** @returns The global mask */
    af::versa<bool, af::c_grid<2> > global_mask() const {
      return global_mask_;
    }

    /** @returns The thresholded value mask */
    af::versa<bool, af::c_grid<2> > value_mask() const {
      return value_mask_;
    }

    /** @returns The final mask of strong pixels. */
    af::versa<bool, af::c_grid<2> > final_mask() const {
      return final_mask_;
    }

  private:
    void init(const af::const_ref<double, af::c_grid<2> > &image,
              const af::const_ref<bool, af::c_grid<2> > &mask,
              const af::const_ref<double, af::c_grid<2> > &gain,
              int2 size,
              double nsig_b,
              double nsig_s,
              double threshold,
              int min_count) {
      // Check the input
      DIALS_ASSERT(threshold >= 0);
      DIALS_ASSERT(nsig_b >= 0 && nsig_s >= 0);
      DIALS_ASSERT(image.accessor().all_eq(mask.accessor()));
      DIALS_ASSERT(image.accessor().all_eq(gain.accessor()));

      // Copy the mask into a temp variable
      af::versa<int, af::c_grid<2> > temp(mask.accessor());
      for (std::size_t i = 0; i < temp.size(); ++i) {
        temp[i] = mask[i] ? 1 : 0;
      }

      // Calculate the masked index_of_dispersion filtered image
      IndexOfDispersionFilterMasked<double> filter(
        image, temp.const_ref(), size, min_count);
      mean_ = filter.mean();
      variance_ = filter.sample_variance();
      cv_ = filter.index_of_dispersion();
      af::versa<int, af::c_grid<2> > count = filter.count();
      temp = filter.mask();

      // Assign pixels to object or background
      cv_mask_ = af::versa<bool, af::c_grid<2> >(image.accessor(), false);
      value_mask_ = af::versa<bool, af::c_grid<2> >(image.accessor(), false);
      final_mask_ = af::versa<bool, af::c_grid<2> >(image.accessor(), false);
      global_mask_ = af::versa<bool, af::c_grid<2> >(image.accessor(), false);
      for (std::size_t i = 0; i < image.size(); ++i) {
        if (temp[i]) {
          double bnd_b = gain[i] + nsig_b * gain[i] * std::sqrt(2.0 / (count[i] - 1));
          double bnd_s = mean_[i] + nsig_s * std::sqrt(gain[i] * mean_[i]);
          cv_mask_[i] = cv_[i] > bnd_b;
          global_mask_[i] = image[i] > threshold;
          value_mask_[i] = image[i] > bnd_s;
          final_mask_[i] = cv_mask_[i] && value_mask_[i] & global_mask_[i];
        }
      }
    }

    af::versa<double, af::c_grid<2> > mean_;
    af::versa<double, af::c_grid<2> > variance_;
    af::versa<double, af::c_grid<2> > cv_;
    af::versa<bool, af::c_grid<2> > global_mask_;
    af::versa<bool, af::c_grid<2> > cv_mask_;
    af::versa<bool, af::c_grid<2> > value_mask_;
    af::versa<bool, af::c_grid<2> > final_mask_;
  };

  /**
   * A class to help debug spot finding by exposing the results of various bits
   * of processing.
   */
  class DispersionExtendedThresholdDebug {
  public:
    /**
     * Do the processing.
     * @param image The image array
     * @param mask The mask array
     * @param size The size of the local window
     * @param nsig_b The background threshold.
     * @param nsig_s The strong pixel threshold
     * @param threshold The global threshold value
     * @param min_count The minimum number of pixels in the local area
     */
    DispersionExtendedThresholdDebug(const af::const_ref<double, af::c_grid<2> > &image,
                                     const af::const_ref<bool, af::c_grid<2> > &mask,
                                     int2 size,
                                     double nsig_b,
                                     double nsig_s,
                                     double threshold,
                                     int min_count) {
      af::versa<double, af::c_grid<2> > gain(image.accessor(), 1.0);
      init(image, mask, gain.const_ref(), size, nsig_b, nsig_s, threshold, min_count);
    }

    /**
     * Do the processing.
     * @param image The image array
     * @param mask The mask array
     * @param size The size of the local window
     * @param nsig_b The background threshold.
     * @param nsig_s The strong pixel threshold
     * @param threshold The global threshold value
     * @param min_count The minimum number of pixels in the local area
     */
    DispersionExtendedThresholdDebug(const af::const_ref<double, af::c_grid<2> > &image,
                                     const af::const_ref<bool, af::c_grid<2> > &mask,
                                     const af::const_ref<double, af::c_grid<2> > &gain,
                                     int2 size,
                                     double nsig_b,
                                     double nsig_s,
                                     double threshold,
                                     int min_count) {
      init(image, mask, gain, size, nsig_b, nsig_s, threshold, min_count);
    }

    /** @returns The mean map */
    af::versa<double, af::c_grid<2> > mean() const {
      return mean_;
    }

    /** @returns The variance map. */
    af::versa<double, af::c_grid<2> > variance() const {
      return variance_;
    }

    /** @returns The index of dispersion map */
    af::versa<double, af::c_grid<2> > index_of_dispersion() const {
      return cv_;
    }

    /** @returns The thresholded index of dispersion mask */
    af::versa<bool, af::c_grid<2> > cv_mask() const {
      return cv_mask_;
    }

    /** @returns The global mask */
    af::versa<bool, af::c_grid<2> > global_mask() const {
      return global_mask_;
    }

    /** @returns The thresholded value mask */
    af::versa<bool, af::c_grid<2> > value_mask() const {
      return value_mask_;
    }

    /** @returns The final mask of strong pixels. */
    af::versa<bool, af::c_grid<2> > final_mask() const {
      return final_mask_;
    }

  private:
    void init(const af::const_ref<double, af::c_grid<2> > &image,
              const af::const_ref<bool, af::c_grid<2> > &mask,
              const af::const_ref<double, af::c_grid<2> > &gain,
              int2 size,
              double nsig_b,
              double nsig_s,
              double threshold,
              int min_count) {
      // Check the input
      DIALS_ASSERT(threshold >= 0);
      DIALS_ASSERT(nsig_b >= 0 && nsig_s >= 0);
      DIALS_ASSERT(image.accessor().all_eq(mask.accessor()));
      DIALS_ASSERT(image.accessor().all_eq(gain.accessor()));

      // Copy the mask into a temp variable
      af::versa<int, af::c_grid<2> > temp(mask.accessor());
      for (std::size_t i = 0; i < temp.size(); ++i) {
        temp[i] = mask[i] ? 1 : 0;
      }

      // Calculate the masked index_of_dispersion filtered image
      IndexOfDispersionFilterMasked<double> filter(
        image, temp.const_ref(), size, min_count);
      mean_ = filter.mean();
      variance_ = filter.sample_variance();
      cv_ = filter.index_of_dispersion();
      af::versa<int, af::c_grid<2> > count = filter.count();
      temp = filter.mask();

      // Assign pixels to object or background
      cv_mask_ = af::versa<bool, af::c_grid<2> >(image.accessor(), false);
      value_mask_ = af::versa<bool, af::c_grid<2> >(image.accessor(), false);
      final_mask_ = af::versa<bool, af::c_grid<2> >(image.accessor(), false);
      global_mask_ = af::versa<bool, af::c_grid<2> >(image.accessor(), false);
      for (std::size_t i = 0; i < image.size(); ++i) {
        if (!temp[i]) {
          /*
           * Masked pixels are set as non-background in order to
           * prevent them from affecting the erosion calculation.
           */
          cv_mask_[i] = true;
        } else {
          double bnd_b = gain[i] + nsig_b * gain[i] * std::sqrt(2.0 / (count[i] - 1));
          cv_mask_[i] = cv_[i] > bnd_b;
          global_mask_[i] = image[i] > threshold;
        }
      }

      // Compute the chebyshev distance to the background (for morphological erosion)
      af::versa<int, af::c_grid<2> > distance(image.accessor(), 0);
      chebyshev_distance(cv_mask_.const_ref(), false, distance.ref());

      // Erode the strong pixel mask and set a mask containing only strong pixels
      std::size_t erosion_distance = std::min(size[0], size[1]);
      af::versa<int, af::c_grid<2> > temp_mask(image.accessor(), 0);
      for (std::size_t i = 0; i < image.size(); ++i) {
        if (temp[i]) {
          value_mask_[i] = distance[i] >= erosion_distance;
          temp_mask[i] = !(cv_mask_[i] && value_mask_[i]);
        }
      }

      // Widen the kernel slightly and compute the mean image without strong pixels
      size[0] += 2;
      size[1] += 2;
      mean2_ = mean_filter_masked(image, temp_mask.ref(), size, 2, false);

      // Compute the final thresholds
      for (std::size_t i = 0; i < image.size(); ++i) {
        if (mask[i]) {
          double bnd_s = mean2_[i] + nsig_s * std::sqrt(gain[i] * mean2_[i]);
          value_mask_[i] = (distance[i] >= erosion_distance) && (image[i] >= bnd_s);
          final_mask_[i] = cv_mask_[i] && value_mask_[i] && global_mask_[i];
        }
        global_mask_[i] = temp_mask[i];
      }
    }

    af::versa<double, af::c_grid<2> > mean_;
    af::versa<double, af::c_grid<2> > mean2_;
    af::versa<double, af::c_grid<2> > variance_;
    af::versa<double, af::c_grid<2> > cv_;
    af::versa<bool, af::c_grid<2> > global_mask_;
    af::versa<bool, af::c_grid<2> > cv_mask_;
    af::versa<bool, af::c_grid<2> > value_mask_;
    af::versa<bool, af::c_grid<2> > final_mask_;
  };

  /**
   * A class to compute the threshold using index of dispersion
   */
  class DispersionExtendedThreshold {
  public:
    /**
     * Enable more efficient memory usage by putting components required for the
     * summed area table closer together in memory
     */
    template <typename T>
    struct Data {
      int m;
      T x;
      T y;
    };

    DispersionExtendedThreshold(int2 image_size,
                                int2 kernel_size,
                                double nsig_b,
                                double nsig_s,
                                double threshold,
                                int min_count)
        : image_size_(image_size),
          kernel_size_(kernel_size),
          nsig_b_(nsig_b),
          nsig_s_(nsig_s),
          threshold_(threshold),
          min_count_(min_count) {
      // Check the input
      DIALS_ASSERT(threshold_ >= 0);
      DIALS_ASSERT(nsig_b >= 0 && nsig_s >= 0);
      DIALS_ASSERT(image_size.all_gt(0));
      DIALS_ASSERT(kernel_size.all_gt(0));

      // Ensure the min counts are valid
      std::size_t num_kernel = (2 * kernel_size[0] + 1) * (2 * kernel_size[1] + 1);
      if (min_count_ <= 0) {
        min_count_ = num_kernel;
      } else {
        DIALS_ASSERT(min_count_ <= num_kernel && min_count_ > 1);
      }

      // Allocate the buffer
      std::size_t element_size = sizeof(Data<double>);
      buffer_.resize(element_size * image_size[0] * image_size[1]);
    }

    /**
     * Compute the summed area tables for the mask, src and src^2.
     * @param src The input array
     * @param mask The mask array
     */
    template <typename T>
    void compute_sat(af::ref<Data<T> > table,
                     const af::const_ref<T, af::c_grid<2> > &src,
                     const af::const_ref<bool, af::c_grid<2> > &mask) {
      // Largest value to consider
      const T BIG = (1 << 24);  // About 16m counts

      // Get the size of the image
      std::size_t ysize = src.accessor()[0];
      std::size_t xsize = src.accessor()[1];

      // Create the summed area table
      for (std::size_t j = 0, k = 0; j < ysize; ++j) {
        int m = 0;
        T x = 0;
        T y = 0;
        for (std::size_t i = 0; i < xsize; ++i, ++k) {
          int mm = (mask[k] && src[k] < BIG) ? 1 : 0;
          m += mm;
          x += mm * src[k];
          y += mm * src[k] * src[k];
          if (j == 0) {
            table[k].m = m;
            table[k].x = x;
            table[k].y = y;
          } else {
            table[k].m = table[k - xsize].m + m;
            table[k].x = table[k - xsize].x + x;
            table[k].y = table[k - xsize].y + y;
          }
        }
      }
    }

    /**
     * Compute the threshold
     * @param src - The input array
     * @param mask - The mask array
     * @param dst The output array
     */
    template <typename T>
    void compute_dispersion_threshold(af::ref<Data<T> > table,
                                      const af::const_ref<T, af::c_grid<2> > &src,
                                      const af::const_ref<bool, af::c_grid<2> > &mask,
                                      af::ref<bool, af::c_grid<2> > dst) {
      // Get the size of the image
      std::size_t ysize = src.accessor()[0];
      std::size_t xsize = src.accessor()[1];

      // The kernel size
      int kxsize = kernel_size_[1];
      int kysize = kernel_size_[0];

      // Calculate the local mean at every point
      for (std::size_t j = 0, k = 0; j < ysize; ++j) {
        for (std::size_t i = 0; i < xsize; ++i, ++k) {
          int i0 = i - kxsize - 1, i1 = i + kxsize;
          int j0 = j - kysize - 1, j1 = j + kysize;
          i1 = i1 < xsize ? i1 : xsize - 1;
          j1 = j1 < ysize ? j1 : ysize - 1;
          int k0 = j0 * xsize;
          int k1 = j1 * xsize;

          // Compute the number of points valid in the local area,
          // the sum of the pixel values and the sum of the squared pixel
          // values.
          double m = 0;
          double x = 0;
          double y = 0;
          if (i0 >= 0 && j0 >= 0) {
            const Data<T> &d00 = table[k0 + i0];
            const Data<T> &d10 = table[k1 + i0];
            const Data<T> &d01 = table[k0 + i1];
            m += d00.m - (d10.m + d01.m);
            x += d00.x - (d10.x + d01.x);
            y += d00.y - (d10.y + d01.y);
          } else if (i0 >= 0) {
            const Data<T> &d10 = table[k1 + i0];
            m -= d10.m;
            x -= d10.x;
            y -= d10.y;
          } else if (j0 >= 0) {
            const Data<T> &d01 = table[k0 + i1];
            m -= d01.m;
            x -= d01.x;
            y -= d01.y;
          }
          const Data<T> &d11 = table[k1 + i1];
          m += d11.m;
          x += d11.x;
          y += d11.y;

          // Compute the thresholds
          dst[k] = false;
          if (!mask[k]) {
            /*
             * Masked pixels are set as non-background in order to
             * prevent them from affecting the erosion calculation.
             */
            dst[k] = true;
          } else if (m >= min_count_ && x >= 0) {
            double a = m * y - x * x - x * (m - 1);
            double c = x * nsig_b_ * std::sqrt(2 * (m - 1));
            dst[k] = (a > c);
          }
        }
      }
    }

    /**
     * Compute the threshold
     * @param src - The input array
     * @param mask - The mask array
     * @param gain - The gain array
     * @param dst The output array
     */
    template <typename T>
    void compute_dispersion_threshold(af::ref<Data<T> > table,
                                      const af::const_ref<T, af::c_grid<2> > &src,
                                      const af::const_ref<bool, af::c_grid<2> > &mask,
                                      const af::const_ref<double, af::c_grid<2> > &gain,
                                      af::ref<bool, af::c_grid<2> > dst) {
      // Get the size of the image
      std::size_t ysize = src.accessor()[0];
      std::size_t xsize = src.accessor()[1];

      // The kernel size
      int kxsize = kernel_size_[1];
      int kysize = kernel_size_[0];

      // Calculate the local mean at every point
      for (std::size_t j = 0, k = 0; j < ysize; ++j) {
        for (std::size_t i = 0; i < xsize; ++i, ++k) {
          int i0 = i - kxsize - 1, i1 = i + kxsize;
          int j0 = j - kysize - 1, j1 = j + kysize;
          i1 = i1 < xsize ? i1 : xsize - 1;
          j1 = j1 < ysize ? j1 : ysize - 1;
          int k0 = j0 * xsize;
          int k1 = j1 * xsize;

          // Compute the number of points valid in the local area,
          // the sum of the pixel values and the num of the squared pixel
          // values.
          double m = 0;
          double x = 0;
          double y = 0;
          if (i0 >= 0 && j0 >= 0) {
            const Data<T> &d00 = table[k0 + i0];
            const Data<T> &d10 = table[k1 + i0];
            const Data<T> &d01 = table[k0 + i1];
            m += d00.m - (d10.m + d01.m);
            x += d00.x - (d10.x + d01.x);
            y += d00.y - (d10.y + d01.y);
          } else if (i0 >= 0) {
            const Data<T> &d10 = table[k1 + i0];
            m -= d10.m;
            x -= d10.x;
            y -= d10.y;
          } else if (j0 >= 0) {
            const Data<T> &d01 = table[k0 + i1];
            m -= d01.m;
            x -= d01.x;
            y -= d01.y;
          }
          const Data<T> &d11 = table[k1 + i1];
          m += d11.m;
          x += d11.x;
          y += d11.y;

          // Compute the thresholds
          dst[k] = false;
          if (!mask[k]) {
            /*
             * Masked pixels are set as non-background in order to
             * prevent them from affecting the erosion calculation.
             */
            dst[k] = true;
          } else if (m >= min_count_ && x >= 0) {
            double a = m * y - x * x;
            double c = gain[k] * x * (m - 1 + nsig_b_ * std::sqrt(2 * (m - 1)));
            dst[k] = (a > c);
          }
        }
      }
    }

    /**
     * Erode the dispersion mask: N.B. this inverts the definition in flight -
     * the purpose of erosion is to extract those pixels which are at least
     * the kernel-width away from the nearest true background pixel.
     * The implementation in this context returns the pixels which are
     * valid for assessing an estimate of the background.
     *
     * @param dst The dispersion mask
     */
    void erode_dispersion_demo(const af::const_ref<bool, af::c_grid<2> > &mask,
                               af::ref<bool, af::c_grid<2> > dst) {
      // array size, slow then fast
      std::size_t ysize = dst.accessor()[0];
      std::size_t xsize = dst.accessor()[1];

      // search distance: N.B. that this is in practice one pixel smaller
      int d = std::min(kernel_size_[0], kernel_size_[1]) - 1;

      // scratch array to store the result, which is then inverted back
      // to the input array
      af::versa<bool, af::c_grid<2> > scr(dst.accessor());

      for (int j = 0, k = 0; j < ysize; ++j) {
        for (std::size_t i = 0; i < xsize; ++i, ++k) {
          // pixel has to be non-background (dst[k]) and valid (mask[k]) to continue
          if (!dst[k] || !mask[k]) {
            scr[k] = false;
            continue;
          }

          // take as a prior that this pixel is non-backgroundy
          bool tmp = true;

          // search over a (2 * d + 1) ** 2 pixel grid for any true background pixel
          for (int _j = -d; _j <= d; _j++) {
            if ((j + _j < 0) || (j + _j >= ysize)) continue;
            for (int _i = -d; _i <= d; _i++) {
              if ((i + _i < 0 || i + _i >= xsize)) continue;
              int _k = (j + _j) * xsize + i + _i;
              if (!dst[_k]) {
                tmp = false;
              }
            }
          }

          scr[k] = tmp;
        }
      }

      // copy mask back, inverting as we go
      for (int j = 0, k = 0; j < ysize; ++j) {
        for (std::size_t i = 0; i < xsize; ++i, ++k) {
          dst[k] = mask[k] && !scr[k];
        }
      }
    }

    /**
     * Erode the dispersion mask: N.B. this inverts the definition in flight
     * @param dst The dispersion mask
     */
    void erode_dispersion_mask(const af::const_ref<bool, af::c_grid<2> > &mask,
                               af::ref<bool, af::c_grid<2> > dst) {
      // The distance array
      af::versa<int, af::c_grid<2> > distance(dst.accessor(), 0);

      // Compute the chebyshev distance to the nearest valid background pixel
      chebyshev_distance(dst, false, distance.ref());

      // The erosion distance
      std::size_t erosion_distance = std::min(kernel_size_[0], kernel_size_[1]);

      // Compute the eroded mask
      for (std::size_t k = 0; k < dst.size(); ++k) {
        if (mask[k]) {
          dst[k] = !(dst[k] && distance[k] >= erosion_distance);
        } else {
          dst[k] = false;
        }
      }
    }

    /**
     * Compute the threshold
     * @param src - The input array
     * @param mask - The mask array
     * @param dst The output array
     */
    template <typename T>
    void compute_final_threshold(af::ref<Data<T> > table,
                                 const af::const_ref<T, af::c_grid<2> > &src,
                                 const af::const_ref<bool, af::c_grid<2> > &mask,
                                 af::ref<bool, af::c_grid<2> > dst) {
      // Get the size of the image
      std::size_t ysize = src.accessor()[0];
      std::size_t xsize = src.accessor()[1];

      // The kernel size
      int kxsize = kernel_size_[1] + 2;
      int kysize = kernel_size_[0] + 2;

      // Calculate the local mean at every point
      for (std::size_t j = 0, k = 0; j < ysize; ++j) {
        for (std::size_t i = 0; i < xsize; ++i, ++k) {
          int i0 = i - kxsize - 1, i1 = i + kxsize;
          int j0 = j - kysize - 1, j1 = j + kysize;
          i1 = i1 < xsize ? i1 : xsize - 1;
          j1 = j1 < ysize ? j1 : ysize - 1;
          int k0 = j0 * xsize;
          int k1 = j1 * xsize;

          // Compute the number of points valid in the local area,
          // the sum of the pixel values and the sum of the squared pixel
          // values.
          double m = 0;
          double x = 0;
          if (i0 >= 0 && j0 >= 0) {
            const Data<T> &d00 = table[k0 + i0];
            const Data<T> &d10 = table[k1 + i0];
            const Data<T> &d01 = table[k0 + i1];
            m += d00.m - (d10.m + d01.m);
            x += d00.x - (d10.x + d01.x);
          } else if (i0 >= 0) {
            const Data<T> &d10 = table[k1 + i0];
            m -= d10.m;
            x -= d10.x;
          } else if (j0 >= 0) {
            const Data<T> &d01 = table[k0 + i1];
            m -= d01.m;
            x -= d01.x;
          }
          const Data<T> &d11 = table[k1 + i1];
          m += d11.m;
          x += d11.x;

          // Compute the thresholds. The pixel is marked True if:
          // 1. The pixel is valid
          // 2. It has 1 or more unmasked neighbours
          // 3. It is within the dispersion masked region
          // 4. It is greater than the global threshold
          // 5. It is greater than the local mean threshold
          //
          // Otherwise it is false
          if (mask[k] && m >= 0 && x >= 0) {
            bool dispersion_mask = !dst[k];
            bool global_mask = src[k] > threshold_;
            double mean = (m >= 2 ? (x / m) : 0);
            bool local_mask = src[k] >= (mean + nsig_s_ * std::sqrt(mean));
            dst[k] = dispersion_mask && global_mask && local_mask;
          } else {
            dst[k] = false;
          }
        }
      }
    }

    /**
     * Compute the threshold
     * @param src - The input array
     * @param mask - The mask array
     * @param dst The output array
     */
    template <typename T>
    void compute_final_threshold(af::ref<Data<T> > table,
                                 const af::const_ref<T, af::c_grid<2> > &src,
                                 const af::const_ref<bool, af::c_grid<2> > &mask,
                                 const af::const_ref<double, af::c_grid<2> > &gain,
                                 af::ref<bool, af::c_grid<2> > dst) {
      // Get the size of the image
      std::size_t ysize = src.accessor()[0];
      std::size_t xsize = src.accessor()[1];

      // The kernel size
      int kxsize = kernel_size_[1] + 2;
      int kysize = kernel_size_[0] + 2;

      // Calculate the local mean at every point
      for (std::size_t j = 0, k = 0; j < ysize; ++j) {
        for (std::size_t i = 0; i < xsize; ++i, ++k) {
          int i0 = i - kxsize - 1, i1 = i + kxsize;
          int j0 = j - kysize - 1, j1 = j + kysize;
          i1 = i1 < xsize ? i1 : xsize - 1;
          j1 = j1 < ysize ? j1 : ysize - 1;
          int k0 = j0 * xsize;
          int k1 = j1 * xsize;

          // Compute the number of points valid in the local area,
          // the sum of the pixel values and the sum of the squared pixel
          // values.
          double m = 0;
          double x = 0;
          if (i0 >= 0 && j0 >= 0) {
            const Data<T> &d00 = table[k0 + i0];
            const Data<T> &d10 = table[k1 + i0];
            const Data<T> &d01 = table[k0 + i1];
            m += d00.m - (d10.m + d01.m);
            x += d00.x - (d10.x + d01.x);
          } else if (i0 >= 0) {
            const Data<T> &d10 = table[k1 + i0];
            m -= d10.m;
            x -= d10.x;
          } else if (j0 >= 0) {
            const Data<T> &d01 = table[k0 + i1];
            m -= d01.m;
            x -= d01.x;
          }
          const Data<T> &d11 = table[k1 + i1];
          m += d11.m;
          x += d11.x;

          // Compute the thresholds. The pixel is marked True if:
          // 1. The pixel is valid
          // 2. It has 1 or more unmasked neighbours
          // 3. It is within the dispersion masked region
          // 4. It is greater than the global threshold
          // 5. It is greater than the local mean threshold
          //
          // Otherwise it is false
          if (mask[k] && m >= 0 && x >= 0) {
            bool dispersion_mask = !dst[k];
            bool global_mask = src[k] > threshold_;
            double mean = (m >= 2 ? (x / m) : 0);
            bool local_mask = src[k] >= (mean + nsig_s_ * std::sqrt(gain[k] * mean));
            dst[k] = dispersion_mask && global_mask && local_mask;
          } else {
            dst[k] = false;
          }
        }
      }
    }

    /**
     * Compute the threshold for the given image and mask.
     * @param src - The input image array.
     * @param mask - The mask array.
     * @param dst - The destination array.
     */
    template <typename T>
    void threshold(const af::const_ref<T, af::c_grid<2> > &src,
                   const af::const_ref<bool, af::c_grid<2> > &mask,
                   af::ref<bool, af::c_grid<2> > dst) {
      // check the input
      DIALS_ASSERT(src.accessor().all_eq(image_size_));
      DIALS_ASSERT(src.accessor().all_eq(mask.accessor()));
      DIALS_ASSERT(src.accessor().all_eq(dst.accessor()));

      // Get the table
      DIALS_ASSERT(sizeof(T) <= sizeof(double));

      // Cast the buffer to the table type
      af::ref<Data<T> > table(reinterpret_cast<Data<T> *>(&buffer_[0]), buffer_.size());

      // compute the summed area table
      compute_sat(table, src, mask);

      // Compute the dispersion threshold. This output is in dst which contains
      // a mask where 1 is valid background and 0 is invalid pixels and stuff
      // above the dispersion threshold
      compute_dispersion_threshold(table, src, mask, dst);

      // Erode the dispersion mask: N.B. this changes in place the definition of
      // dst from "pixels that are not background" to "pixels that are background"
      erode_dispersion_mask(mask, dst);

      // Compute the summed area table again now excluding the threshold pixels
      // (which are set to false in dst)
      compute_sat(table, src, dst);

      // Compute the final threshold
      compute_final_threshold(table, src, mask, dst);
    }

    /**
     * Compute the threshold for the given image and mask.
     * @param src - The input image array.
     * @param mask - The mask array.
     * @param gain - The gain array
     * @param dst - The destination array.
     */
    template <typename T>
    void threshold_w_gain(const af::const_ref<T, af::c_grid<2> > &src,
                          const af::const_ref<bool, af::c_grid<2> > &mask,
                          const af::const_ref<double, af::c_grid<2> > &gain,
                          af::ref<bool, af::c_grid<2> > dst) {
      // check the input
      DIALS_ASSERT(src.accessor().all_eq(image_size_));
      DIALS_ASSERT(src.accessor().all_eq(mask.accessor()));
      DIALS_ASSERT(src.accessor().all_eq(gain.accessor()));
      DIALS_ASSERT(src.accessor().all_eq(dst.accessor()));

      // Get the table
      DIALS_ASSERT(sizeof(T) <= sizeof(double));

      // Cast the buffer to the table type
      af::ref<Data<T> > table((Data<T> *)&buffer_[0], buffer_.size());

      // compute the summed area table
      compute_sat(table, src, mask);

      // Compute the dispersion threshold. This output is in dst which contains
      // a mask where 1 is valid background and 0 is invalid pixels and stuff
      // above the dispersion threshold
      compute_dispersion_threshold(table, src, mask, gain, dst);

      // Erode the dispersion mask
      erode_dispersion_mask(mask, dst);

      // Compute the summed area table again now excluding the threshold pixels
      compute_sat(table, src, dst);

      // Compute the final threshold
      compute_final_threshold(table, src, mask, gain, dst);
    }

  private:
    int2 image_size_;
    int2 kernel_size_;
    double nsig_b_;
    double nsig_s_;
    double threshold_;
    int min_count_;
    std::vector<char> buffer_;
  };

}}  // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_IMAGE_THRESHOLD_LOCAL_H */
