# Import the necessary libraries
import argparse
import tkinter as tk
from tkinter import filedialog

import numpy as np
from PIL import Image, ImageTk
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import subprocess

from matplotlib.figure import Figure

from config import config
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

parser = argparse.ArgumentParser()
parser.add_argument('--createdb', action='store_true', default=False, help='Whether the database should be created')

args = parser.parse_args()

# Parse config file
paths, opt = config('config')

if args.createdb:
    # Run makeblastdb command to create a database for BLAST
    command = f"{paths['BLAST']}/makeblastdb -in {paths['DATABASE_PATH']} -dbtype nucl"
    subprocess.run(command, shell=True)

# Specify the path to your local database
DATABASE_PATH = paths['DATABASE_PATH']


# Function to run BLAST on a given sequence against the local database
def blast_sequence(sequence):
    blastn_cline = NcbiblastnCommandline(cmd=paths['BLAST'] + '/blastn', query=sequence, db=DATABASE_PATH, outfmt=5,
                                         out="blast_results.xml")
    stdout, stderr = blastn_cline()
    result_handle = open("blast_results.xml")
    blast_record = NCBIXML.read(result_handle)
    return blast_record


min_hsp = (float('inf'), None)


# Function to execute the BLAST run and display results
def run_blast():
    global min_hsp
    # Get the sequence from the entry box
    # sequence = seq_entry.get()
    filename = filedialog.askopenfilename(initialdir="./",
                                          title="Select a File",
                                          filetypes=(("Seq files",
                                                      "*.seq*"),
                                                     ("fasta files",
                                                      "*.fa*"),
                                                     ("all files",
                                                      "*.*")))
    # Run BLAST with the sequence
    blast_record = blast_sequence(filename)
    # Display results
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect <= min_hsp[0]:
                min_hsp = (hsp.expect, hsp, alignment)
            results_text.insert(tk.END, f"****Alignment****\n")
            results_text.insert(tk.END, f"sequence: {alignment.title}\n")
            results_text.insert(tk.END, f"length: {alignment.length}\n")
            results_text.insert(tk.END, f"e value: {hsp.expect}\n")
            results_text.insert(tk.END, hsp.query[0:75] + "...\n")
            results_text.insert(tk.END, hsp.match[0:75] + "...\n")
            results_text.insert(tk.END, hsp.sbjct[0:75] + "...\n")
            titles_text.insert(tk.END, f"sequence: {alignment.title}\n")


def _delta(x, y):
    return 0 if x == y else 1


def _M(seq1, seq2, i, j, k):
    return sum(_delta(x, y) for x, y in zip(seq1[i:i + k], seq2[j:j + k]))


def _makeMatrix(seq1, seq2, k):
    n = len(seq1)
    m = len(seq2)
    return [[_M(seq1, seq2, i, j, k) for j in range(m - k + 1)] for i in range(n - k + 1)]


def plot():
    # the figure that will contain the plot
    fig = Figure(figsize=(5, 5),
                 dpi=100)
    hsp = min_hsp[1]

    print(min_hsp)
    win = win_entry.get()
    if win:
        win = int(win)
    else:
        win = 10
    if hsp:
        query = hsp.query
        sbjct = hsp.sbjct

        # list of squares
        y = np.array(_makeMatrix(query, sbjct, win))
        print(len(query), len(sbjct))
        # adding the subplot
        plot1 = fig.add_subplot(111)
        fig.suptitle(min_hsp[2].title)

        # plotting the graph
        print(y)
        plot1.imshow(y)

    # creating the Tkinter canvas
    # containing the Matplotlib figure
    canvas = FigureCanvasTkAgg(fig,
                               master=Tick)
    canvas.draw()

    # placing the canvas on the Tkinter window
    canvas.get_tk_widget().grid(column=1, row=0)


# Create a Tkinter window
Tick = tk.Tk()

# Load an image and display it
image1 = Image.open("TICK_logo.png")
test = ImageTk.PhotoImage(image1)
label1 = tk.Label(image=test)
label1.image = test
label1.grid(row=0, column=1)

# Create a frame within the window
mainframe = tk.Frame(Tick, bg='black', padx=10, pady=10)
mainframe.grid(column=0, row=0, sticky=(tk.N, tk.W, tk.E, tk.S))

# Create label and entry box for sequence input
# seq_label = tk.Label(mainframe, text="Enter sequence:", bg='white')
# seq_label.grid(column=0, row=1, sticky=tk.W)
# seq_entry = tk.Entry(mainframe, bg="white", width=50)
# seq_entry.grid(column=1, row=1, sticky=tk.E)

# Create "Run BLAST" button
run_button = tk.Button(mainframe, text="Run BLAST", command=run_blast)
run_button.grid(column=0, row=3, columnspan=3)

# Create "Plot" button
run_button = tk.Button(mainframe, text="Dotplot", command=plot)
run_button.grid(column=0, row=9, columnspan=1)
win_label = tk.Label(mainframe, text="Enter dotplot window (default: 10):", bg='white')
win_label.grid(column=1, row=9, sticky=tk.W)
win_entry = tk.Entry(mainframe, bg="white", width=10)
win_entry.grid(column=2, row=9, sticky=tk.E)


# Create a text box for displaying BLAST results
results_text = tk.Text(mainframe, height=10, width=100, bg="gray")
results_text.grid(column=0, row=5, columnspan=2)

# Create a text box for displaying BLAST results titles
titles_text = tk.Text(mainframe, height=10, width=100, bg="gray")
titles_text.grid(column=0, row=7, columnspan=2)


# Create a scrollbar for the text box
scrollbar = tk.Scrollbar(mainframe, command=results_text.yview)
scrollbar.grid(row=5, column=2, sticky='nsew')
results_text['yscrollcommand'] = scrollbar.set

# Start the Tkinter event loop
Tick.mainloop()
