# -*- coding: utf-8 -*-
#
# Copyright (c) 2023 Alberto Casagrande <alberto.casagrande@uniud.it>
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
            0: RACES.CellEventType.DEATH,
            1: RACES.CellEventType.DUPLICATION,
            2: RACES.CellEventType.EPIGENETIC_SWITCH,
            3: RACES.CellEventType.GENOTYPE_MUTATION
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


class TestGenotype(unittest.TestCase):
    def test_init(self):
        try:
            RACES.Genotype("A", [RACES.EpigeneticRates(0.01, 0.01)])
        except BaseException:
            self.fail('RACES.Genotype("A", '
                      + '[RACES.EpigeneticRates(0.01, 0.01)]) '
                      + ' raised an unexpected exeception!')

        try:
            RACES.Genotype("A", [[0.01, 0.01]])
        except BaseException:
            self.fail('RACES.Genotype("A", [[0.01, 0.01]]) '
                      + ' raised an unexpected exeception!')

        try:
            RACES.Genotype("A", [])
        except BaseException:
            self.fail('RACES.Genotype("A", []) '
                      + ' raised an unexpected exeception!')

        try:
            RACES.Genotype("A", [[0.01, 0.01], [0.01, 0.01]])
        except BaseException:
            self.fail('RACES.Genotype("A", [[0.01, 0.01], '
                      + '[0.01, 0.01]]) raised an unexpected exeception!')

        with self.assertRaises(Exception):
            RACES.Genotype("A", 2)

        with self.assertRaises(Exception):
            RACES.Genotype("A", ['a'])

        with self.assertRaises(Exception):
            RACES.Genotype("A", [[0.01, 0.01, 2]])

    def test_properties(self):
        A = RACES.Genotype("A", [[0.01, 0.01]])

        self.assertEqual(A.name, "A")
        self.assertEqual(A.num_of_promoters, 1)

        B = RACES.Genotype("B", [])

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
        A = RACES.Genotype("A", [[0.01, 0.01]])
        try:
            A.set_rates("-", {RACES.CellEventType.DEATH: 0.1,
                              RACES.CellEventType.DUPLICATION: 0.2})
        except BaseException:
            self.fail('A.set_rates("-", {RACES.CellEventType.DEATH: 0.1, '
                      + 'RACES.CellEventType.DUPLICATION: 0.2})'
                      + ' raised an unexpected exeception!')

        with self.assertRaises(Exception):
            A.set_rates("", {RACES.CellEventType.DEATH: 0.1,
                             RACES.CellEventType.DUPLICATION: 0.2})

        B = RACES.Genotype("B", [])

        try:
            B.set_rates("", {RACES.CellEventType.DEATH: 0.1,
                             RACES.CellEventType.DUPLICATION: 0.2})
        except BaseException:
            self.fail('B.set_rates("", {RACES.CellEventType.DEATH: 0.1, '
                      + 'RACES.CellEventType.DUPLICATION: 0.2})'
                      + ' raised an unexpected exeception!')

        with self.assertRaises(Exception):
            B.set_rates("-", {RACES.CellEventType.DEATH: 0.1,
                              RACES.CellEventType.DUPLICATION: 0.2})

    def test_get_rate(self):
        A = RACES.Genotype("A", [[0.01, 0.01]])
        A.set_rates("-", {RACES.CellEventType.DEATH: 0.1,
                          RACES.CellEventType.DUPLICATION: 0.2})
        try:
            self.assertEqual(A.get_rate("-", RACES.CellEventType.DEATH), 0.1)
            self.assertEqual(A.get_rate("-", RACES.CellEventType.DUPLICATION),
                             0.2)
            self.assertEqual(A.get_rate("-",
                                        RACES.CellEventType.EPIGENETIC_SWITCH),
                             0.0)
        except BaseException:
            self.fail('RACES.Genotype.set_rate() '
                      + ' raised an unexpected exeception!')

        with self.assertRaises(Exception):
            self.assertEqual(A.get_rate("", RACES.CellEventType.DEATH), 0.01)

        B = RACES.Genotype("B", [])
        B.set_rates("", {RACES.CellEventType.DEATH: 0.1,
                         RACES.CellEventType.DUPLICATION: 0.2})
        try:
            self.assertEqual(B.get_rate("", RACES.CellEventType.DEATH), 0.1)
            self.assertEqual(B.get_rate("", RACES.CellEventType.DUPLICATION),
                             0.2)
            self.assertEqual(B.get_rate("",
                                        RACES.CellEventType.EPIGENETIC_SWITCH),
                             0.0)
        except BaseException:
            self.fail('RACES.Genotype.set_rate() '
                      + ' raised an unexpected exeception!')

        with self.assertRaises(Exception):
            self.assertEqual(B.get_rate("-", RACES.CellEventType.DEATH), 0.01)


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

    def test_add_genotype(self):
        sim = RACES.Simulation()

        A = RACES.Genotype("A", [[0.01, 0.01]])
        A.set_rates("-", {RACES.CellEventType.DEATH: 0.1,
                          RACES.CellEventType.DUPLICATION: 0.2})
        A.set_rates("+", {RACES.CellEventType.DEATH: 0.01,
                          RACES.CellEventType.DUPLICATION: 0.02})

        sim.set_tissue("Liver", [100, 100])

        try:
            # now the simulation has a tissue
            sim.add_genotype(A)

        except BaseException:
            self.fail('sim.add_genotype(A) raised an unexpected exeception!')

        with self.assertRaises(Exception):
            sim.add_genotype("A")

    def test_place_cell(self):
        sim = RACES.Simulation()

        A = RACES.Genotype("A", [[0.01, 0.01]])
        A.set_rates("-", {RACES.CellEventType.DEATH: 0.1,
                          RACES.CellEventType.DUPLICATION: 0.2})
        A.set_rates("+", {RACES.CellEventType.DEATH: 0.01,
                          RACES.CellEventType.DUPLICATION: 0.02})

        sim.set_tissue("Liver", [100, 100])

        sim.add_genotype(A)

        try:
            sim.place_cell(A, "-", [50, 50])

        except BaseException:
            self.fail('sim.place_cell(A, "-", [50, 50]) raised ' +
                      'an unexpected exeception!')

        with self.assertRaises(Exception):
            sim.place_cell(A, "-", [50, 150])

    def test_schedule_genotype_mutation(self):
        sim = RACES.Simulation()

        sim.set_tissue("Liver", [100, 100])

        A = RACES.Genotype("A", [[0.01, 0.01]])
        A.set_rates("-", {RACES.CellEventType.DEATH: 0.1,
                          RACES.CellEventType.DUPLICATION: 0.2})
        A.set_rates("+", {RACES.CellEventType.DEATH: 0.01,
                          RACES.CellEventType.DUPLICATION: 0.02})
        sim.add_genotype(A)

        B = RACES.Genotype("B", [[0.01, 0.01]])
        B.set_rates("-", {RACES.CellEventType.DEATH: 0.1,
                          RACES.CellEventType.DUPLICATION: 0.3})
        B.set_rates("+", {RACES.CellEventType.DEATH: 0.1,
                          RACES.CellEventType.DUPLICATION: 0.45})
        sim.add_genotype(B)

        C = RACES.Genotype("C", [])
        C.set_rates("", {RACES.CellEventType.DEATH: 0.1,
                         RACES.CellEventType.DUPLICATION: 0.3})
        sim.add_genotype(C)

        try:
            sim.schedule_genotype_mutation(A, B, 70)
        except BaseException:
            self.fail('sim.schedule_genotype_mutation(A, B, 70) raised'
                      + ' an unexpected exeception!')

        with self.assertRaises(Exception):
            # methylation signature incompatible
            sim.schedule_genotype_mutation(A, C, 70)

        try:
            # methylation signature compatible
            sim.schedule_genotype_mutation(C, B, 70)
        except BaseException:
            self.fail('sim.schedule_genotype_mutation(C, B, 70) raised'
                      + ' an unexpected exeception!')

    def test_run_up_to(self):
        sim = RACES.Simulation()

        sim.set_tissue("Liver", [100, 100])

        A = RACES.Genotype("A", [[0.01, 0.01]])
        A.set_rates("-", {RACES.CellEventType.DEATH: 0.1,
                          RACES.CellEventType.DUPLICATION: 0.2})
        A.set_rates("+", {RACES.CellEventType.DEATH: 0.01,
                          RACES.CellEventType.DUPLICATION: 0.02})
        sim.add_genotype(A)

        B = RACES.Genotype("B", [[0.01, 0.01]])
        B.set_rates("-", {RACES.CellEventType.DEATH: 0.1,
                          RACES.CellEventType.DUPLICATION: 0.3})
        B.set_rates("+", {RACES.CellEventType.DEATH: 0.1,
                          RACES.CellEventType.DUPLICATION: 0.45})
        sim.add_genotype(B)
        sim.schedule_genotype_mutation(A, B, 70)

        sim.storage_enabled = False

        sim.run_up_to(0, quiet=True)

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
