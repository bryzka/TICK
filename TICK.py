# Import the necessary libraries
import argparse
import tkinter as tk
from tkinter import filedialog, messagebox, simpledialog, ttk
from typing import List
import logging
import Bio
import time
import threading
import numpy as np
from Bio import SeqIO
from PIL import Image, ImageTk
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML,NCBIWWW
import subprocess
from matplotlib.figure import Figure
from config import config
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

logger = logging.getLogger()
logger.setLevel(logging.WARN)
console_handler = logging.StreamHandler()
logger.addHandler(console_handler)

parser = argparse.ArgumentParser()
parser.add_argument('config', action='store', help='Path to the config file')
parser.add_argument('--createdb', '-c', action='store_true', default=False,
                    help='Whether the database should be created')
parser.add_argument('--web', '-w', action='store_true', default=False, help='should the web blast be used')

args = parser.parse_args()

# Parse config file
paths, opt = config(args.config)

if args.createdb:
    # Run makeblastdb command to create a database for BLAST
    command = f"{paths['BLAST']}/makeblastdb -in {paths['DATABASE_PATH']} -dbtype nucl"
    subprocess.run(command, shell=True)

# Specify the path to your local database
DATABASE_PATH = paths['DATABASE_PATH']


def get_orgs(path: str) -> List[str]:
    """
    Gets lits of organisms to blast sequences against from a file
    :param path: path to file containing list of organisms
    :return: list of organisms
    """
    with open(path, 'r') as f:
        ret = f.readlines()
    return ret


# Function to run BLAST on a given sequence against the local database
def blast_sequence(sequence, start, end):
    if start != 0 or end != 1:
        sequ = SeqIO.read(sequence, "fasta")
        logging.info(sequ)
        sequ.seq = sequ.seq[start: -end]
        SeqIO.write(sequ, sequence+'_trimmed', 'fasta')
    blastn_cline = NcbiblastnCommandline(cmd=paths['BLAST'] + '/blastn', query=sequence, db=DATABASE_PATH, outfmt=5,
                                         out="blast_results.xml")
    stdout, stderr = blastn_cline()
    logging.error(stderr)
    result_handle = open("blast_results.xml")
    blast_record = NCBIXML.parse(result_handle)
    return blast_record


# Function to execute BLAST online
def blast_online(fasta, start, end):
    orgs = get_orgs(paths['ORGS'])
    fasta = SeqIO.read(fasta, "fasta").seq[start: -end]
    logging.info(fasta)
    blast_ret = NCBIWWW.qblast('blastn', 'nt', fasta, entrez_query="[organism] OR ".join(orgs))
    blast_ret = NCBIXML.parse(blast_ret)
    return blast_ret


min_hsp = (float('inf'), None)

def progress_bar():
    progress = ttk.Progressbar(mainframe, orient="horizontal", length=300, mode="indeterminate")
    progress.grid(column=0, row=4, columnspan=3)
    progress.start()
    time.sleep(5)
    progress.stop()

# Function to execute the BLAST run locally and display results
def run_blast():
    global min_hsp

    filename = filedialog.askopenfilename(initialdir="./",
                                          title="Select sequence which you want to blast against database",
                                          filetypes=(
                                              ("fasta files",
                                               "*.fa*"),
                                              ("Seq files",
                                               "*.seq*"),
                                              ("all files",
                                               "*.*"),
                                                     ))

    user_trim = messagebox.askyesno("Question", "Do you want to trim this sequence?")
    start, end = 0, 1
    if user_trim:
        start = simpledialog.askinteger("Input", "How much to trim from the start:", parent=Tick)
        end = simpledialog.askinteger("Input", "How much to trim from the end:", parent=Tick)

    progress_thread = threading.Thread(target=progress_bar)
    progress_thread.start()

    # Run BLAST with the sequence
    try:
        if paths['BLAST'] not in {'WWW', ''} and not args.web:
            blast_record = blast_sequence(filename, start, end)
        else:
            blast_record = blast_online(filename, start, end)
    except Bio.Application.ApplicationError:
        raise Bio.Application.ApplicationError('Error')
    # Display results
    logging.info(blast_record)
    for record in blast_record:
        for alignment in record.alignments:
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
    else:
        results_text.insert(tk.END, 'No more alignments')
        titles_text.insert(tk.END, 'No more alignments')


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

    logging.info(min_hsp)
    win = win_entry.get()
    if win:
        win = int(win)
    else:
        win = 10
    if hsp:
        query = hsp.query
        sbjct = hsp.sbjct

        # creating matrix
        y = np.array(_makeMatrix(query, sbjct, win))
        logging.info(len(query), len(sbjct))
        # adding the subplot
        plot1 = fig.add_subplot(111)
        fig.suptitle(min_hsp[2].title)

        # plotting the graph
        logging.info(y)
        plot1.imshow(y)

    # creating the Tkinter canvas
    # containing the Matplotlib figure
    canvas = FigureCanvasTkAgg(fig,
                               master=Tick)
    canvas.draw()

    # placing the canvas on the Tkinter window
    canvas.get_tk_widget().grid(column=1, row=0)


# ... pozostały kod ...

# Create a Tkinter window
Tick = tk.Tk()
Tick.title("TICK")  
Tick.geometry("1200x700") 
Tick.configure(bg="lightgray") 

icon_image = ImageTk.PhotoImage(file="TICK_logo.jpg") 
Tick.iconphoto(False, icon_image)

# Load an image and display it
image1 = Image.open("TICK_logo.jpg").resize((200, 200))  
test = ImageTk.PhotoImage(image1)
label1 = tk.Label(image=test, bg="lightgray")
label1.grid(row=0, column=2, padx=20, pady=20)

# Custom styles for widgets
style = ttk.Style()
style.configure("TButton", font=("Arial", 12), padding=10)
style.configure("TLabel", font=("Arial", 12))
style.configure("TEntry", font=("Arial", 12))

# Create a frame within the window
mainframe = ttk.Frame(Tick, padding="20 20 20 20")
mainframe.grid(column=0, row=0, sticky=(tk.N, tk.W, tk.E, tk.S))

# Create "Run BLAST" button
run_button = ttk.Button(mainframe, text="Select sequence to run BLAST", command=run_blast)
run_button.grid(column=0, row=3, columnspan=3, pady=10)

# Create "Plot" button
plot_button = ttk.Button(mainframe, text="Dotplot", command=plot)
plot_button.grid(column=0, row=9, columnspan=1, pady=10)

# Window label and entry for dotplot
win_label = ttk.Label(mainframe, text="Enter dotplot window (default: 10):")
win_label.grid(column=1, row=9, sticky=tk.E, padx=10)  
win_entry = ttk.Entry(mainframe, width=10, font=("Arial", 12))
win_entry.grid(column=2, row=9, sticky=tk.W, padx=10)  


# Create a text box for displaying BLAST results
results_text = tk.Text(mainframe, height=10, width=100, bg="white", padx=10, pady=10, borderwidth=2, relief="groove")
results_text.grid(column=0, row=5, columnspan=3, pady=10)

# Create a text box for displaying BLAST results titles
titles_text = tk.Text(mainframe, height=10, width=100, bg="white", padx=10, pady=10, borderwidth=2, relief="groove")
titles_text.grid(column=0, row=7, columnspan=3, pady=10)

# Create a scrollbar for the text box
scrollbar = ttk.Scrollbar(mainframe, command=results_text.yview)
scrollbar.grid(row=5, column=3, sticky='nsew')
results_text['yscrollcommand'] = scrollbar.set

# Start the Tkinter event loop
Tick.mainloop()
