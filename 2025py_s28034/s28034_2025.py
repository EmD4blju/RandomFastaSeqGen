# @ Mikołaj Warda
# The program aims to generate random DNA sequences and save them to a .fasta file with provided user settings.

import random
import logging
from pathlib import Path

# ORIGINAL
# Lack of a Class
# MODIFIED: Overall the code I think is made cleaner and easier to follow. 
# Code is now a few steps closer for a universal use.

class FastaSeqProcessor:
    """Manages operation for the whole random DNA sequence generation process with statistics calculation.
    """
    def __init__(self):
        """During init of FastaSeqProcessor it will ask about details such as: DNA sequence length, username, sequence ID and sequence description.
        """
        self.setup_logger() # sets up logger 
        self.length = None
        while self.length is None or self.length <= 0: # Asks user for a sequence length until he provides a valid one.
            self.length = int(input('Pass sequence length: '))
            if self.length <= 0:
                print(f'Length must be additive. Passed: {self.length}')
        self.seq_id = input("Pass sequence ID: ").strip() # Asks user for a sequence id.
        self.description = input("Pass sequence description: ").strip() # Asks user for a sequence description.
        self.name = input("Pass username: ").strip() # Asks user for a username.
        self.sequence = self.generate_dna_sequence() # Generates a random ATCG DNA sequence.
        self.logger.info(msg=f'FastaSeqProcessor created\tInput values: {self.length, self.seq_id, self.description, self.name}') # Displays information about created instance

    def setup_logger(self):
        """Sets up logger for readability.
        """
        self.logger = logging.getLogger(name=__name__) # Creates logger with the module's name which is __main__ in this scenario
        self.logger.setLevel(level=logging.INFO) # Sets the logger's level to be infomational
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s') # Sets display format for the logger
        stream_handler = logging.StreamHandler() # Creates a handler to write logging records to a stream.
        stream_handler.setFormatter(fmt=formatter) # Sets the formatter for that handler.
        self.logger.addHandler(hdlr=stream_handler) # Adds handler for the logger.
        
# ORIGINAL
# def generate_dna_sequence(length):
#     return ''.join(random.choices(['A', 'C', 'G', 'T'], k=length))
# MODIFIED: is now a method of an instance and is being invoked in initializer to work seamlessly at FastaSeqProcessor creation. 
    def generate_dna_sequence(self) -> str:
        """Generates random DNA sequence of a specified length.
        Returns:
            str: Generated DNA sequence. 
        """
        return ''.join(random.choices(['A', 'C', 'G', 'T'], k=self.length)) # Joins a randomly generated list of ACGT together and returns it.

# ORIGINAL
# def insert_name_into_sequence(sequence, name):
#     position = random.randint(0, len(sequence))
#     return sequence[:position] + name + sequence[position:]
# MODIFIED: is now a method of an instance and does not change the sequence. 
# Think of the method as an attribute which is being calculated when needed.
# Therefore calculation of statistics is simpler.
# It is needed only when saving to a file.
    def get_sequence_with_name_inserted(self) -> str:
        """Returns sequence with randomly inserted specified username.
        Returns:
            str: Modified DNA sequence with name inserted.
        """
        position = random.randint(0, len(self.sequence)) # Gets a random integer for username insertion into DNA sequence.
        return self.sequence[:position] + self.name + self.sequence[position:] # Returns a DNA sequence with username inserted beginning at position index.

# ORIGINAL
# def calculate_statistics(sequence):
#     pure_sequence = ''.join([n for n in sequence if n in 'ACGT'])
#     length = len(pure_sequence)
#     counts = {nucleotide: pure_sequence.count(nucleotide) for nucleotide in 'ACGT'}
#     percentages = {n: (counts[n] / length * 100) for n in 'ACGT'}
#     cg = counts['C'] + counts['G']
#     at = counts['A'] + counts['T']
#     cg_at_ratio = (cg / at * 100) if at != 0 else 0
#     return percentages, cg_at_ratio
# MODIFIED: is now a method of an instance and is being invoked when we want to display statistics of a sequence in a display_stats() method.
# Is made simpler due to the held DNA sequence in memory as an instance's attribute.
# The sequence is not influeced by the inserted name anymore.
    def get_stats(self) -> tuple[dict, float]:
        """Calculates statistics for the specified sequence.
        Returns:
            tuple[dict, float]: A tuple containing (dictionary of ACGT percentages and CG ratio).
        """
        counts = {nucleotide: self.sequence.count(nucleotide) for nucleotide in 'ACGT'} # Creates a dictionary/map to hold counts for each nucleotide in the sequence.
        percentages = {n: (counts[n] / self.length * 100) for n in 'ACGT'} # Creates a dictionary/map to hold percentage values of nucleotides occurances in the sequence.
        cg = counts['C'] + counts['G'] # Sum of C counts and G counts.
        at = counts['A'] + counts['T'] # Sum of A counts and T counts.
        cg_at_ratio = (cg / at * 100) if at != 0 else 0 # Calculates ratio of CG to AT. To prevent division by zero the condition is implemented.
        return percentages, cg_at_ratio # Returns calculated values.

    def display_stats(self) -> None:
        """Displays stats of generated DNA sequence.
        """
        percentages, cg_at_ratio = self.get_stats() # Calculates statistics for the DNA sequence.
        self.logger.info(msg="Sequence statistics:") # Displays information about statistics.
        for n in 'ACGT':
            self.logger.info(msg=f"{n}: {percentages[n]:.1f}%")
        self.logger.info(msg=f"%CG: {percentages['C'] + percentages['G']:.1f}")
        self.logger.info(msg=f"CG to AT ratio: {cg_at_ratio:.2f}")

# ORIGINAL
# def save_to_fasta(filename, header, sequence):
#     with open(filename, 'w') as file:
#         file.write(f">{header}\n")
#         file.write(sequence + "\n")
# MODIFIED: is now a method of an instance. Does not require any arguments passed because those are given by the user at the initialization.
# Saves .fasta files in a dedicated saved_fasta_sequences directory.
# Pathlib is used to enable work within most OS's (MacOS, Windows, GNU).
    def save_to_fasta(self) -> None:
        """Saves a sequence to a .fasta file in saved_fasta_sequences directory
        """
        file_path = Path('saved_fasta_sequences', f'{self.seq_id}.fasta') # Creates a file path.
        with open(file_path, 'w') as file: # Try with resources to automatically close the stream.
            file.write(f">{self.seq_id} {self.description}\n") # Writes the header in .fasta format.
            file.write(self.get_sequence_with_name_inserted() + "\n") # Writes the DNA sequence.
        self.logger.info(msg=f"\nSequence was saved to file: {file_path}") # Displays information about successfull save.
    
    def __str__(self):
        return f'Sequence: {self.sequence},\nSequence id: {self.seq_id},\nUsername: {self.name},\nDescription: {self.description}' # Just overriden __str__ (like Java toString()).

# ORIGINAL
# def main():
#     try:
#         length = int(input("Podaj długość sekwencji: "))
#         if length <= 0:
#             raise ValueError("Długość musi być liczbą dodatnią.")
#     except ValueError as e:
#         print(f"Błąd: {e}")
#         return

#     seq_id = input("Podaj ID sekwencji: ").strip()
#     description = input("Podaj opis sekwencji: ").strip()
#     name = input("Podaj imię: ").strip()

#     dna_sequence = generate_dna_sequence(length)
#     dna_with_name = insert_name_into_sequence(dna_sequence, name)
#     percentages, cg_at_ratio = calculate_statistics(dna_with_name)

#     filename = f"{seq_id}.fasta"
#     header = f"{seq_id} {description}"
#     save_to_fasta(filename, header, dna_with_name)

#     print(f"\nSekwencja została zapisana do pliku {filename}")
#     print("Statystyki sekwencji:")
#     for n in 'ACGT':
#         print(f"{n}: {percentages[n]:.1f}%")
#     print(f"%CG: {percentages['C'] + percentages['G']:.1f}")
#     print(f"Stosunek CG do AT: {cg_at_ratio:.2f}")
# MODIFIED: main became smaller due to major modifications made with the rest of the code and process happening within FastaSeqProcessor framework.
def main():
    fsp = FastaSeqProcessor() # Initialization of FastaSeqProcessor instance.
    fsp.display_stats() # Statistic display.
    fsp.save_to_fasta() # Saving to .fasta file.

if __name__ == "__main__": # Runs the code.
    main()