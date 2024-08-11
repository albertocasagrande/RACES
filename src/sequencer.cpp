/**
 * @file sequencer.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements sequencer models
 * @version 1.0
 * @date 2024-08-11
 *
 * @copyright Copyright (c) 2023-2024
 *
 * MIT License
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


#include "sequencer.hpp"

namespace RACES
{

namespace Sequencers
{

ConstantQualityScoreModel::ConstantQualityScoreModel():
    quality_score(RACES::Sequencers::BasicSequencer::max_phred_code)
{}

ConstantQualityScoreModel::ConstantQualityScoreModel(const uint8_t quality_score):
    quality_score(std::min(quality_score, BasicSequencer::max_phred_code))
{}

namespace Illumina
{

double BasicQualityScoreModel::default_mean(const size_t& posisiton)
{
    const double library_const = static_cast<double>(100)/150;
    constexpr double mean_a1 = 0.0845729; //!< The mean linear coefficient before the peak
    const double mean_a2 = -0.0146604*library_const; //!< The mean linear coefficient after the peak
    constexpr double mean_b = 35.1498051;  //!< The mean quality score in first position
    
    constexpr uint16_t mean_peak_pos = 7; //!< The position of the mean peak

    if (posisiton<mean_peak_pos) {
        return mean_a1 * posisiton + mean_b;
    }

    return mean_a2*posisiton + (mean_peak_pos*(mean_a1-mean_a2)+mean_b);
}

double BasicQualityScoreModel::default_stddev(const size_t& posisiton)
{
    const double library_const = static_cast<double>(100)/150;

    const double stddev_quad_coeff = 9.661484e-06 * library_const * library_const;   //!< The standard deviation quadratic coefficient
    constexpr double stddev_base_coeff = 1.731231e-01;   //!< The standard deviation zero-degree coefficient

    return stddev_quad_coeff*posisiton*posisiton + stddev_base_coeff;
}

BasicQualityScoreModel::BasicQualityScoreModel(const BasicQualityScoreModel::EstimatingFunctionType& mean_qual_function,
                                               const BasicQualityScoreModel::EstimatingFunctionType& stddev_qual_function):
    mean(mean_qual_function), stddev(stddev_qual_function)
{}

}   // Illumina

}   // Sequencers

}   // RACES
