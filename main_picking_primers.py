# PICK PRIMERS

class Primer:

    def __init__(self, dna_string):
        """
        Initializes the Primer object with the provided DNA sequence.

        Parameters:
        dna_string (str): The DNA sequence for the primer.
        """
        self.seq = dna_string
        self.len = len(dna_string)
        self.gc = self.gc_content()
        self.tm = self.primer_tm()
        self.min17 = self.min_large_17()
        self.end = dna_string[-1]
        self.comp = self.primer_comp()

    def gc_content(self):
        """
        Calculates and returns the GC content of the primer.

        Returns:
        float: GC content as a percentage of the total length of the primer.
        """
        return 100 * sum(1 for nuc in self.seq if nuc in 'GC') / self.len

    def min_large_17(self):
        """
        Checks if the length of the primer is at least 17 base pairs.

        Returns:
        bool: True if the primer is 17 bp or longer, False otherwise.
        """
        return self.len >= 17

    def primer_tm(self):
        """
        Calculates and returns the melting temperature (Tm) of the primer.

        Returns:
        int: The Tm of the primer calculated using the formula 4*(G+C) + 2*(A+T).
        """
        return 4 * (self.seq.count('C') + self.seq.count('G')) + 2 * (self.seq.count('A') + self.seq.count('T'))

    def primer_comp(self):
        """
        Generates and returns the complementary DNA sequence of the primer.

        This method takes the DNA sequence of the primer stored in `self.seq` and
        generates its complementary sequence. The complementary bases are determined
        based on the following pairs: A-T, T-A, C-G, G-C.

        Returns:
        str: The complementary DNA sequence of the primer.
        """
        complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join([complement_dict[i] for i in self.seq])[::-1]

    def no_hairpin(self):
        """
        Checks if the primer can form a hairpin structure.

        Parameters:
        min_length (int): Minimum length of the palindromic sequence to consider.

        Returns:
        bool: True if no hairpin can form, False otherwise.
        """
        complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        list_of_hairpin_trials_1 =[]
        list_of_hairpin_trials_2 = []

        for i in range(self.len - 3):
            list_of_hairpin_trials_1.append(self.seq[i:i+4])
            list_of_hairpin_trials_2.append(self.seq[i:i + 4])

        #list_of_hairpin_trials_2.reverse()

        for i in range(len(list_of_hairpin_trials_2)):
            list_of_hairpin_trials_2.insert(i, ''.join([complement_dict[nuc] for nuc in list_of_hairpin_trials_2[i]]))

        return (list_of_hairpin_trials_1, list_of_hairpin_trials_2)


    @staticmethod
    def no_primer_dimer(primer_1, primer_2):
        """
        Checks if the two primers form primer-dimers.

        Parameters:
        primer_1 (Primer): The first primer object.
        primer_2 (Primer): The second primer object.

        Returns:
        bool: True if the two primers do not form primer-dimers, False otherwise.
        """
        complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        last_nuc1, last_nuc2 = primer_1.seq[-4:], primer_2.seq.comp[-4:]
        number_of_comp = 0
        for nuc1, nuc2 in zip(last_nuc1, last_nuc2):
            if nuc1 == nuc2:
                number_of_comp += 1
        return number_of_comp <= 2

    @staticmethod
    def gc_content_test(primer_1, primer_2):
        """
        Checks if the GC content of both primers is at least 50%.

        Parameters:
        primer_1 (Primer): The first primer object.
        primer_2 (Primer): The second primer object.

        Returns:
        bool: True if both primers have at least 50% GC content, False otherwise.
        """
        return primer_1.gc >= 50 and primer_2.gc >= 50

    @staticmethod
    def anneling_temperature_check(primer_1, primer_2):
        """
        Checks if the difference in melting temperature (Tm) between two primers is within 5°C.

        Parameters:
        primer_1 (Primer): The first primer object.
        primer_2 (Primer): The second primer object.

        Returns:
        bool: True if the Tm difference is less than or equal to 5°C, False otherwise.
        """
        return abs(primer_1.tm - primer_2.tm) <= 5

    @staticmethod
    def is_end_g_or_c(primer_1, primer_2):
        """
        Checks if both primers end with a G or C nucleotide.

        Parameters:
        primer_1 (Primer): The first primer object.
        primer_2 (Primer): The second primer object.

        Returns:
        bool: True if both primers end with G or C, False otherwise.
        """
        return primer_1.end in ('G', 'C') and primer_2.end in ('G', 'C')
