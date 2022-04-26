NAME
INTERCAAT (Interface Contact definition with Adaptable Atom Types software)

SYNOPSIS
This program uses a PDB file to identify the residues present in the interface 
between a query chain and an interacting chain(s)

PRE-REQUISITES
1. Python 3.8 or newer is required.
2. Numpy is required
3. Qhull is recommended since it is much faster. Qhull can be downloaded from here: http://www.qhull.org/download/.
4. If you wish not to use Qhull, the Python package SciPy is needed.

INTERCAAT was developed using numpy 1.21.1 and scipy 1.7.0.
They can be installed in $HOME with the below commands:
pip install numpy --user   
pip install scipy --user 

TEST
We included a sample PDB and the expected results from running the test commands below:
python intercaat.py -pdb 1cph.pdb -qc B -ic A 
python intercaat.py -pdb 1cph.pdb -qc A,6,7,8,9 -ic B -di no -cc no -sr 1.7 

CONFIGURATION
For ease of use the program is set to run with Python implementation as a default.
To run with qhull, you must set "run_python_version = no" and specify the "qvoronoi_path" and
"qhull_path" path to run the qvoronoi executable in the intercaat_config.ini file.

HELP
The help function for INTERCAAT can be displayed by typing either intercaat.py or intercaat.py -h into 
your terminal. After doing so, instructions for how to run this program along with an 
example for each switch will be shown.

OUTPUT
The main output of this program displays every atomic interaction between the query chain 
and the interacting chain(s), the distance between the interacting atoms, and the assigned 
atom classes. The compatibility matrix, if displayed, shows each interface residue in the 
query chain and its corresponding quantity of atomic interactions.

REFERENCE
INTERCAAT: identifying interface residues be-tween macromolecules

AUTHORS
Steven Grudman	steven.grudman@einsteinmed.org
Eduardo Fajardo	eduardo@fiserlab.org
Andras Fiser    andras.fiser@einsteinmed.org

