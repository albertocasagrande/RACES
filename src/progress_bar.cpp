/**
 * @file progress_bar.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Implement a progress bar
 * @version 0.5
 * @date 2023-08-05
 * 
 * @copyright Copyright (c) 2023
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

#include "progress_bar.hpp"

namespace Races
{

namespace UI
{

ProgressBar::ProgressBar():
    last_update(std::chrono::system_clock::now()), percentage(0), message(), updated(false)
#ifdef WITH_INDICATORS
    , indicator()
#endif // WITH_INDICATORS
{
    using namespace std::chrono_literals;

    update_interval = 1s;

#ifdef WITH_INDICATORS
    using namespace indicators;

    indicator = new indicators::ProgressBar{
        option::BarWidth{40},
        option::Start{" ["},
        option::Fill{"█"},
        option::Lead{"█"},
        option::Remainder{"-"},
        option::End{"]"},
        option::ForegroundColor{Color::yellow},
        option::ShowElapsedTime{true},
        option::ShowPercentage{true},
        option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
    };
#endif  // WITH_INDICATORS

}

ProgressBar& ProgressBar::set_message(const std::string& message)
{

    this->message = message;

#ifdef WITH_INDICATORS
    using namespace indicators;

    if (percentage != 100 || !updated) {
        indicator->set_option(option::PostfixText{message});
        indicator->set_progress(percentage);
    }
#endif  // WITH_INDICATORS

    updated = true;

    return *this;
}

ProgressBar& ProgressBar::set_progress(const unsigned int percentage, const std::string& message)
{
    using namespace std::chrono;

    this->percentage = std::min(percentage, static_cast<unsigned int>(100));
    this->message = message;

    updated = false;

    const auto curr_time = system_clock::now();
    if (duration_cast<seconds>(curr_time-last_update).count()>=1) {
        last_update = curr_time;

        set_message(message);
    }

    return *this;
}


ProgressBar& ProgressBar::set_progress(const unsigned int percentage)
{
    using namespace std::chrono;

    this->percentage = percentage;

    updated = false;

    const auto curr_time = system_clock::now();
    if (duration_cast<seconds>(curr_time-last_update).count()>=1) {
        last_update = curr_time;

        updated = true;

#ifdef WITH_INDICATORS
        indicator->set_progress(percentage);
#endif  // WITH_INDICATORS
    }

    return *this;
}

const unsigned int& ProgressBar::get_progress() const
{
    return percentage;
}

void ProgressBar::show_console_cursor()
{
#ifdef WITH_INDICATORS
    indicators::show_console_cursor(true);
#endif  // WITH_INDICATORS
}

void ProgressBar::hide_console_cursor()
{
#ifdef WITH_INDICATORS
    indicators::show_console_cursor(false);
#endif  // WITH_INDICATORS
}

ProgressBar::~ProgressBar()
{
#ifdef WITH_INDICATORS
    set_message(this->message);

    if (percentage != 100) {
        indicator->mark_as_completed();
    }

    show_console_cursor();

    delete indicator;
#endif  // WITH_INDICATORS
}

}  // UI

} // Races