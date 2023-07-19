# -*- coding: utf-8 -*-
#
# Copyright (c) 2023 Alberto Casagrande <acasagrande@units.it>
#
# MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import unittest

import RACES


class TestCellEventType(unittest.TestCase):

    def test_values(self):
        self.assertEqual(RACES.CellEventType.values, {
            0: RACES.CellEventType.DIE,
            1: RACES.CellEventType.DUPLICATE,
            2: RACES.CellEventType.EPIGENETIC_EVENT,
            3: RACES.CellEventType.DUPLICATION_AND_EPIGENETIC_EVENT,
            4: RACES.CellEventType.PASSENGER_MUTATION,
            5: RACES.CellEventType.DRIVER_SOMATIC_MUTATION
        })


class TestEpigeneticRates(unittest.TestCase):
    def test_init(self):
        try:
            RACES.EpigeneticRates(0.2, 0.1)
        except BaseException:
            self.fail("RACES.EpigeneticRates(0.2, 0.1) unexpectedly"
                      + " raised an exeception!")

        with self.assertRaises(Exception):
            RACES.EpigeneticRates()

        with self.assertRaises(Exception):
            RACES.EpigeneticRates(-0.2, 0.1)

        with self.assertRaises(Exception):
            RACES.EpigeneticRates(0.2, -0.1)

        with self.assertRaises(Exception):
            RACES.EpigeneticRates(1.2, 0.1)

        with self.assertRaises(Exception):
            RACES.EpigeneticRates(0.2, 1.1)

    def test_get(self):
        e_rates = RACES.EpigeneticRates(0.2, 0.1)

        self.assertEqual(e_rates.get_methylation_rate(), 0.2)
        self.assertEqual(e_rates.get_demethylation_rate(), 0.1)

    def test_set(self):
        e_rates = RACES.EpigeneticRates(0.2, 0.1)

        e_rates.set_methylation_rate(0.3)
        e_rates.set_demethylation_rate(0.4)

        self.assertEqual(e_rates.get_methylation_rate(), 0.3)
        self.assertEqual(e_rates.get_demethylation_rate(), 0.4)

        with self.assertRaises(Exception):
            e_rates.set_methylation_rate(1.3)

        with self.assertRaises(Exception):
            e_rates.set_demethylation_rate(1.3)

        with self.assertRaises(Exception):
            e_rates.set_methylation_rate(-0.3)

        with self.assertRaises(Exception):
            e_rates.set_demethylation_rate(-0.3)


class TestSomaticGenotype(unittest.TestCase):
    def test_init(self):
        try:
            RACES.SomaticGenotype("A", [RACES.EpigeneticRates(0.01, 0.01)])
        except BaseException:
            self.fail('RACES.SomaticGenotype("A", '
                      + '[RACES.EpigeneticRates(0.01, 0.01)]) '
                      + ' raised an unexpected exeception!')

        try:
            RACES.SomaticGenotype("A", [[0.01, 0.01]])
        except BaseException:
            self.fail('RACES.SomaticGenotype("A", [[0.01, 0.01]]) '
                      + ' raised an unexpected exeception!')

        try:
            RACES.SomaticGenotype("A", [])
        except BaseException:
            self.fail('RACES.SomaticGenotype("A", []) '
                      + ' raised an unexpected exeception!')

        try:
            RACES.SomaticGenotype("A", [[0.01, 0.01], [0.01, 0.01]])
        except BaseException:
            self.fail('RACES.SomaticGenotype("A", [[0.01, 0.01], '
                      + '[0.01, 0.01]]) raised an unexpected exeception!')

        with self.assertRaises(Exception):
            RACES.SomaticGenotype("A", 2)

        with self.assertRaises(Exception):
            RACES.SomaticGenotype("A", ['a'])

        with self.assertRaises(Exception):
            RACES.SomaticGenotype("A", [[0.01, 0.01, 2]])

    def test_properties(self):
        A = RACES.SomaticGenotype("A", [[0.01, 0.01]])

        self.assertEqual(A.name, "A")
        self.assertEqual(A.num_of_promoters, 1)

        B = RACES.SomaticGenotype("B", [])

        self.assertEqual(B.name, "B")
        self.assertEqual(B.id, A.id+1)
        self.assertEqual(B.num_of_promoters, 0)

        with self.assertRaises(Exception):
            A.name = "test"

        with self.assertRaises(Exception):
            A.id = 3

        with self.assertRaises(Exception):
            A.num_of_promoters = 7

    def test_set_rates(self):
        A = RACES.SomaticGenotype("A", [[0.01, 0.01]])
        try:
            A.set_rates("-", {RACES.CellEventType.DIE: 0.1,
                              RACES.CellEventType.DUPLICATE: 0.2})
        except BaseException:
            self.fail('A.set_rates("-", {RACES.CellEventType.DIE: 0.1, '
                      + 'RACES.CellEventType.DUPLICATE: 0.2})'
                      + ' raised an unexpected exeception!')

        with self.assertRaises(Exception):
            A.set_rates("", {RACES.CellEventType.DIE: 0.1,
                             RACES.CellEventType.DUPLICATE: 0.2})

        B = RACES.SomaticGenotype("B", [])

        try:
            B.set_rates("", {RACES.CellEventType.DIE: 0.1,
                             RACES.CellEventType.DUPLICATE: 0.2})
        except BaseException:
            self.fail('B.set_rates("", {RACES.CellEventType.DIE: 0.1, '
                      + 'RACES.CellEventType.DUPLICATE: 0.2})'
                      + ' raised an unexpected exeception!')

        with self.assertRaises(Exception):
            B.set_rates("-", {RACES.CellEventType.DIE: 0.1,
                              RACES.CellEventType.DUPLICATE: 0.2})

    def test_get_rate(self):
        A = RACES.SomaticGenotype("A", [[0.01, 0.01]])
        A.set_rates("-", {RACES.CellEventType.DIE: 0.1,
                          RACES.CellEventType.DUPLICATE: 0.2})
        try:
            self.assertEqual(A.get_rate("-", RACES.CellEventType.DIE), 0.1)
            self.assertEqual(A.get_rate("-", RACES.CellEventType.DUPLICATE),
                             0.2)
            self.assertEqual(A.get_rate("-",
                                        RACES.CellEventType.EPIGENETIC_EVENT),
                             0.0)
        except BaseException:
            self.fail('RACES.SomaticGenotype.set_rate() '
                      + ' raised an unexpected exeception!')

        with self.assertRaises(Exception):
            self.assertEqual(A.get_rate("", RACES.CellEventType.DIE), 0.01)

        B = RACES.SomaticGenotype("B", [])
        B.set_rates("", {RACES.CellEventType.DIE: 0.1,
                         RACES.CellEventType.DUPLICATE: 0.2})
        try:
            self.assertEqual(B.get_rate("", RACES.CellEventType.DIE), 0.1)
            self.assertEqual(B.get_rate("", RACES.CellEventType.DUPLICATE),
                             0.2)
            self.assertEqual(B.get_rate("",
                                        RACES.CellEventType.EPIGENETIC_EVENT),
                             0.0)
        except BaseException:
            self.fail('RACES.SomaticGenotype.set_rate() '
                      + ' raised an unexpected exeception!')

        with self.assertRaises(Exception):
            self.assertEqual(B.get_rate("-", RACES.CellEventType.DIE), 0.01)


class TestSimulation(unittest.TestCase):
    def test_init(self):
        try:
            RACES.Simulation()
        except BaseException:
            self.fail('RACES.Simulation() raised an unexpected exeception!')

        try:
            RACES.Simulation(5)
        except BaseException:
            self.fail('RACES.Simulation(5) raised an unexpected exeception!')

        try:
            RACES.Simulation(5, 1)
        except BaseException:
            self.fail('RACES.Simulation(5, 1) raised an unexpected'
                      + ' exeception!')

    def test_set_tissue(self):
        sim = RACES.Simulation()

        try:
            sim.set_tissue("Liver", [100, 100])
        except BaseException:
            self.fail('sim.set_tissue("Liver", [100, 100]) raised'
                      + ' an unexpected exeception!')

        try:
            sim.set_tissue("Liver", [100, 100, 100])
        except BaseException:
            self.fail('sim.set_tissue("Liver", [100, 100, 100]) raised'
                      + ' an unexpected exeception!')

        with self.assertRaises(Exception):
            sim.set_tissue("Liver", [100])

        with self.assertRaises(Exception):
            sim.set_tissue("Liver", [])

        with self.assertRaises(Exception):
            sim.set_tissue("Liver", [100, 1, 1, 1])

        with self.assertRaises(Exception):
            sim.set_tissue("Liver", [100, 100, -100])

        with self.assertRaises(Exception):
            sim.set_tissue("Liver", [100, -100])

    def test_get_time(self):
        sim = RACES.Simulation()

        self.assertEqual(sim.get_time(), 0)

    def test_add_species(self):
        sim = RACES.Simulation()

        A = RACES.SomaticGenotype("A", [[0.01, 0.01]])
        A.set_rates("-", {RACES.CellEventType.DIE: 0.1,
                          RACES.CellEventType.DUPLICATE: 0.2})
        A.set_rates("+", {RACES.CellEventType.DIE: 0.01,
                          RACES.CellEventType.DUPLICATE: 0.02})

        with self.assertRaises(Exception):
            # no tissue has been associated to the simulation yet
            sim.add_species(A)

        sim.set_tissue("Liver", [100, 100])

        try:
            # now the simulation has a tissue
            sim.add_species(A)

        except BaseException:
            self.fail('sim.add_species(A) raised an unexpected exeception!')

        with self.assertRaises(Exception):
            sim.add_species("A")

    def test_add_species(self):
        sim = RACES.Simulation()

        A = RACES.SomaticGenotype("A", [[0.01, 0.01]])
        A.set_rates("-", {RACES.CellEventType.DIE: 0.1,
                          RACES.CellEventType.DUPLICATE: 0.2})
        A.set_rates("+", {RACES.CellEventType.DIE: 0.01,
                          RACES.CellEventType.DUPLICATE: 0.02})

        sim.set_tissue("Liver", [100, 100])

        sim.add_species(A)

        try:
            sim.add_cell(A, "-", [50, 50])

        except BaseException:
            self.fail('sim.add_species(A) raised an unexpected exeception!')

        with self.assertRaises(Exception):
            sim.add_cell(A, "-", [50, 150])

    def test_add_somatic_mutation(self):
        sim = RACES.Simulation()

        sim.set_tissue("Liver", [100, 100])

        A = RACES.SomaticGenotype("A", [[0.01, 0.01]])
        A.set_rates("-", {RACES.CellEventType.DIE: 0.1,
                          RACES.CellEventType.DUPLICATE: 0.2})
        A.set_rates("+", {RACES.CellEventType.DIE: 0.01,
                          RACES.CellEventType.DUPLICATE: 0.02})
        sim.add_species(A)

        B = RACES.SomaticGenotype("B", [[0.01, 0.01]])
        B.set_rates("-", {RACES.CellEventType.DIE: 0.1,
                          RACES.CellEventType.DUPLICATE: 0.3})
        B.set_rates("+", {RACES.CellEventType.DIE: 0.1,
                          RACES.CellEventType.DUPLICATE: 0.45})
        sim.add_species(B)

        C = RACES.SomaticGenotype("C", [])
        C.set_rates("", {RACES.CellEventType.DIE: 0.1,
                         RACES.CellEventType.DUPLICATE: 0.3})
        sim.add_species(C)

        try:
            sim.add_somatic_mutation(A, B, 70)
        except BaseException:
            self.fail('sim.add_somatic_mutation(A, B, 70) raised'
                      + ' an unexpected exeception!')

        with self.assertRaises(Exception):
            # methylation signature incompatible
            sim.add_somatic_mutation(A, C, 70)

        try:
            # methylation signature compatible
            sim.add_somatic_mutation(C, B, 70)
        except BaseException:
            self.fail('sim.add_somatic_mutation(C, B, 70) raised'
                      + ' an unexpected exeception!')

    def test_run_up_to(self):
        sim = RACES.Simulation()

        sim.set_tissue("Liver", [100, 100])

        A = RACES.SomaticGenotype("A", [[0.01, 0.01]])
        A.set_rates("-", {RACES.CellEventType.DIE: 0.1,
                          RACES.CellEventType.DUPLICATE: 0.2})
        A.set_rates("+", {RACES.CellEventType.DIE: 0.01,
                          RACES.CellEventType.DUPLICATE: 0.02})
        sim.add_species(A)

        B = RACES.SomaticGenotype("B", [[0.01, 0.01]])
        B.set_rates("-", {RACES.CellEventType.DIE: 0.1,
                          RACES.CellEventType.DUPLICATE: 0.3})
        B.set_rates("+", {RACES.CellEventType.DIE: 0.1,
                          RACES.CellEventType.DUPLICATE: 0.45})
        sim.add_species(B)
        sim.add_somatic_mutation(A, B, 70)

        sim.run_up_to(100, logging=False, quiet=True)

    def test_death_activation_level(self):
        sim = RACES.Simulation()

        try:
            value = sim.death_activation_level
        except BaseException:
            self.fail('reading sim.death_activation_level raised'
                      + ' an unexpected exeception!')

        self.assertTrue(isinstance(value, int))

        try:
            sim.death_activation_level = 100
        except BaseException:
            self.fail('setting sim.death_activation_level raised'
                      + ' an unexpected exeception!')

        self.assertEqual(sim.death_activation_level, 100)


if __name__ == '__main__':
    unittest.main()
