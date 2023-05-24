# Import the necessary libraries
import tkinter as tk
from PIL import Image, ImageTk
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import subprocess

# Run makeblastdb command to create a database for BLAST
command = "makeblastdb -in local_database.fasta -dbtype nucl"
subprocess.run(command, shell=True)

# Specify the path to your local database
DATABASE_PATH = "local_database.fasta" 

# Function to run BLAST on a given sequence against the local database
def blast_sequence(sequence):
    blastn_cline = NcbiblastnCommandline(query=sequence, db=DATABASE_PATH, outfmt=5, out="blast_results.xml")
    stdout, stderr = blastn_cline()
    result_handle = open("blast_results.xml")
    blast_record = NCBIXML.read(result_handle)
    return blast_record

# Function to execute the BLAST run and display results
def run_blast():
    # Get the sequence from the entry box
    sequence = seq_entry.get()
    # Run BLAST with the sequence
    blast_record = blast_sequence(sequence)
    # Display results
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            results_text.insert(tk.END, f"****Alignment****\n")
            results_text.insert(tk.END, f"sequence: {alignment.title}\n")
            results_text.insert(tk.END, f"length: {alignment.length}\n")
            results_text.insert(tk.END, f"e value: {hsp.expect}\n")
            results_text.insert(tk.END, hsp.query[0:75] + "...\n")
            results_text.insert(tk.END, hsp.match[0:75] + "...\n")
            results_text.insert(tk.END, hsp.sbjct[0:75] + "...\n")

# Create a Tkinter window
Tick = tk.Tk()

# Load an image and display it
image1 = Image.open("TICK_logo.jpg")
test = ImageTk.PhotoImage(image1)
label1 = tk.Label(image=test)
label1.image = test
label1.grid(row=0, column=1)

# Create a frame within the window
mainframe = tk.Frame(Tick, bg='black', padx=10, pady=10)
mainframe.grid(column=0, row=0, sticky=(tk.N, tk.W, tk.E, tk.S))

# Create label and entry box for sequence input
seq_label = tk.Label(mainframe, text="Enter sequence:", bg='white')
seq_label.grid(column=0, row=1, sticky=tk.W)
seq_entry = tk.Entry(mainframe, bg="white", width=50)
seq_entry.grid(column=1, row=1, sticky=tk.E)

# Create "Run BLAST" button
run_button = tk.Button(mainframe, text="Run BLAST", command=run_blast)
run_button.grid(column=0, row=3, columnspan=3)

# Create a text box for displaying BLAST results
results_text = tk.Text(mainframe, height=10, width=50, bg="gray")
results_text.grid(column=0, row=5, columnspan=2)

# Create a scrollbar for the text box
scrollbar = tk.Scrollbar(mainframe, command=results_text.yview)
scrollbar.grid(row=5, column=2, sticky='nsew')
results_text['yscrollcommand'] = scrollbar.set

# Start the Tkinter event loop
Tick.mainloop()
