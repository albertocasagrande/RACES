/**
 * @file simulator.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Main file for the simulator
 * @version 0.1
 * @date 2023-05-30
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

#include <iostream> 

#include <vector>

#include "logger.hpp"
#include "simulator.hpp"

#include "SDL_plot.hpp"

int main()
{
    using namespace Races;

    std::vector<DriverGenotype> genotypes;

    genotypes.push_back(DriverGenotype("A",{
            {CellEventType::DIE, 0.002},{CellEventType::DUPLICATE, 0.1}}));
    genotypes.push_back(DriverGenotype("A",{
            {CellEventType::DIE, 0.001},{CellEventType::DUPLICATE, 0.3}},true));

    genotypes.push_back(DriverGenotype("B",{
            {CellEventType::DIE, 0.002},{CellEventType::DUPLICATE, 0.2}}));
    genotypes.push_back(DriverGenotype("B",{
            {CellEventType::DIE, 0.001},{CellEventType::DUPLICATE, 0.02}},true));

    Tissue tissue("Liver", 1000,1000);

    for (const auto& gen: genotypes) {
        tissue.add_species(gen);
    }

    tissue.add_driver_epigenetic_mutation(genotypes[0].get_id(),genotypes[1].get_id(), 0.001);
    tissue.add_driver_epigenetic_mutation(genotypes[1].get_id(),genotypes[0].get_id(), 0.001);
    tissue.add_driver_epigenetic_mutation(genotypes[2].get_id(),genotypes[3].get_id(), 0.001);
    tissue.add_driver_epigenetic_mutation(genotypes[3].get_id(),genotypes[2].get_id(), 0.001);

    tissue.add_driver_somatic_mutation(genotypes[0].get_id(),genotypes[2].get_id(), 100);
    tissue.add_driver_somatic_mutation(genotypes[1].get_id(),genotypes[3].get_id(), 100);

    tissue.add(genotypes[0].get_id(), {250, 500}, 0);
    tissue.add(genotypes[2].get_id(), {750, 500}, 0);

    BasicSimulator<BasicLogger, UI::SDLWindow> simulator(tissue);
    //BasicSimulator<BasicLogger> simulator(tissue);

    simulator.snapshot_interval = 150000;

    simulator.run_up_to(300000);

    return 0;
}
