# This program is used to identify a reasonably-optimized set of genomes for a certain taxa to include in an MS-GF+ database, such that all peptides in the samples for that taxa will be identified.
# This program reads in a set of MS-GF+ search results that were run on a database consisting only of sequences from the taxa of interest.
# The program then identifies a set of species that would hit all of the signficant peptide hits in a subsequent MS-GF+ search.
# Finally, the program writes out the list of species to stdout, and makes a new fasta file with the sequences for all the genomes that should be included in a final database.
# ARGUMENT1 - Path to the directory with TSV search results
# ARGUMENT2 - Fasta file used to run the initial MS-GF+ search
# ARGUMENT3 - The full path, including the name of the fasta file, where the program should output results
# ex: python DatabaseSpeciesOptimizer.py Results/CommunityV3 Sequences/CompiledDatabases/GardnerellaOnly.fasta DatabaseOptimization/gardnerella_optimization.csv
import sys, csv
from pathlib import Path

# OBJECT AND FUNCTION DEFINITIONS
# Object tracking a specific peptide sequence, the names of genomes currently covering it, and how many times the specific peptide came up in MS-GF+ search results.
class Peptide:
	def __init__(self, sequence):
		self.seq = sequence
		self.coveredBy = set()
		self.timesFound = 1

# Object representing the genome of a species that tracks what sample peptides it covers.
class Genome:
	def __init__(self, name):
		self.name = name
		self.coveredPeptides = set()
	
	@property
	def numPeptides(self):
		return len(self.coveredPeptides)

# Object used to store information about a protein
class Protein:
	def __init__(self, entry):
		self.id = self.extractProteinID(entry)
		self.taxa = self.extractSpecies(entry)
		self.name = self.extractProteinName(entry)
		self.data = entry[entry.find('OX='):entry.find('\n')]
		self.sequence = entry[entry.find('\n') + 1:]
		self.hitPeptides = set()
	
	# This function takes a single protein entry from a .fasta file and returns all the species for that entry
	def extractSpecies(self, entry):
		stringStart = entry.find('OS=') + 3
		stringEnd = entry.find(' OX=')
		return entry[stringStart:stringEnd].split(';')

	# This function takes a single protein entry from a .fasta file and returns the protein identifier for that entry.
	def extractProteinID(self, entry):
		stringEnd = entry.find(' ')
		return entry[0:stringEnd]

	# This function takes a single protein entry from a .fasta file and returns the protein identifier for that entry.
	def extractProteinName(self, entry):
		stringStart = entry.find(' ') + 1
		stringEnd = entry.find(' OS=')
		return entry[stringStart:stringEnd]
	
	# Get the fasta entry for this protein, including all taxa that encode it
	def getEntry(self):
		taxaEntry = self.taxa[0]
		if len(self.taxa) > 1:
			for i in range(1, len(self.taxa)):
				taxaEntry = f'{taxaEntry};{self.taxa[i]}'
		return f'>{self.id} {self.name} OS={taxaEntry} {self.data}\n{self.sequence}\n'
	
	def equals(self, otherProtein):
		return self.id == otherProtein.id

# Object used for tying protein data to an organism.
# Uses a hash set as a data backend so that finding a protein by ID is O(1).
class ProtRef:
	# Create a new ProtRef object using the specified dbFile Path. If disallowRepeats = True, the constructor will throw an exception if there are proteins with identical IDs.
	def __init__(self, dbFilePath, disallowRepeats=True):
		self.db = {}
		rawText = dbFilePath.read_text()
		dataList = rawText.split('\n>')
		dataList[0] = dataList[0][1:]
		del rawText
		for sequence in dataList:
			if sequence.find('Contaminant_') != -1:
				continue
			newProt = Protein(sequence)
			if disallowRepeats and newProt.id in self.db.keys():
				raise Exception('There are at least two proteins in the database with ID "' + newProt.id + '". Protein IDs must be unique.')
			self.db[newProt.id] = newProt
	
	# Retrieves a protein from the reference database by its ID.
	def getProt(self, id):
		try:
			return self.db[id]
		except:
			raise Exception('No protein with ID "' + id + '" was found in the database.')
	
	# Get the name of the taxa for the protein with the specified ID.
	def getTaxaFor(self, id):
		try:
			return self.db[id].taxa
		except:
			raise Exception('Could not find the taxa for protein "' + id + '" because no such protein exists in the database.')
	
	# Returns a randomly ordered list of all the taxa with proteins in this reference database
	def getTaxaList(self):
		taxaSet = set()
		for x in self.db.values():
			for y in x.taxa:
				taxaSet.add(y)
		return list(taxaSet)

# Returns a list of protein hits for the spectrum
def getProteinHitList(row):
	proteins = row[9].split(';')
	toReturn = []
	for i in range(0, len(proteins)):
		endIndex = proteins[i].find('(')
		toReturn.append(proteins[i][0:endIndex])
	return toReturn

# Takes a protein ID and returns what type of organism the protein belongs to.
# Can return the strings 'human', 'decoy', 'contaminant', 'first' (when the "protein" is actually the first row of the TSV file), or 'nonhuman'
def determineProteinType(ids):
	for protID in ids:
		if protID.find('XXX') != -1:
			return 'decoy'
		elif protID.find('Contaminant_') != -1:
			return 'contaminant'
		elif protID.find('HUMAN') != -1:
			return 'human'
	return 'nonhuman'

# MAIN PROGRAM THREAD
# Check user inputs
if len(sys.argv) < 2:
	raise Exception('You must provide a directory with TSV search results')
tsvPath = Path(sys.argv[1])
if not tsvPath.exists():
	raise Exception(f'Could not find directory {sys.argv[1]}')
if len(sys.argv) < 3:
	raise Exception('You must provide a .fasta database file')
dbPath = Path(sys.argv[2])
if not dbPath.exists():
	raise Exception(f'Could not find the specified .fasta file {sys.argv[2]}')
if not dbPath.suffix == '.fasta':
	raise Exception('Database file must be in fasta format')
if len(sys.argv) < 4:
	raise Exception('You must enter the name of a fasta file where output database can be written')
outPath = Path(sys.argv[3])
if not outPath.suffix == '.fasta':
	raise Exception('Output file must be in fasta format')

# Make reference objects
print('Reading database...')
ref = ProtRef(dbPath)
speciesList = {}
for taxa in ref.getTaxaList():
	speciesList[taxa] = Genome(taxa)

# Get a list of all peptides in the samples, associate them with genomes
print('Collecting peptide hits...')
peptides = {}
specCount = 0
for resFile in [x for x in tsvPath.iterdir() if x.suffix == '.tsv']:
	with resFile.open(mode='r') as infile:
		tsvReader = csv.reader(infile, delimiter='\t')
		for row in tsvReader:
			if row[14] == 'QValue':
				continue
			if float(row[14]) > 0.01:
				break
			hitList = getProteinHitList(row)
			if determineProteinType(hitList) != 'nonhuman':
				continue
			sequence = row[8]
			if sequence in peptides.keys():
				peptides[sequence].timesFound += 1
				specCount += 1
			else:
				peptides[sequence] = Peptide(sequence)
				specCount += 1
			for hit in hitList:
				ref.getProt(hit).hitPeptides.add(sequence)
				for species in ref.getProt(hit).taxa:
					peptides[sequence].coveredBy.add(species)
					speciesList[species].coveredPeptides.add(sequence)
print(f'Total spec count: {specCount}')

# Only include peptides whose proteins had more than one unique peptide identified across all samples
peptidesToInclude = set()
elCount1 = 0
for prot in ref.db.values():
	if len(prot.hitPeptides) > 1:
		for pep in prot.hitPeptides:
			peptidesToInclude.add(pep)
for pep in [pep for pep in peptidesToInclude if len(peptides[pep].coveredBy) == 1]:
	elCount1 += 1
print(f'Total number of unique peptides: {str(len(peptides))}')
print(f'Number of unique peptides after filtering: {str(len(peptidesToInclude))}')
print(f'Number of filtered peptides covered by only 1 genome: {str(elCount1)}')

# Make the list of genomes that will hit every included peptide
print('Determining optimal database taxa list...')
genomeList = []
prevCovered = set()
workingList = []
for genome in speciesList.values():
	workingList.append(genome)
bestDiff = 1000000000000000000000
while bestDiff > 0:
	bestGenome = -1
	for j in range(len(workingList)):
		tempCovered = prevCovered.copy()
		tempCovered.update(workingList[j].coveredPeptides)
		delta = len(peptidesToInclude.difference(tempCovered))
		if delta < bestDiff:
			bestDiff = delta
			bestGenome = j
	prevCovered.update(workingList[bestGenome].coveredPeptides)
	genomeList.append(workingList[bestGenome])
	del workingList[bestGenome]

# Write the final list of taxa to stdout
print(f'Building database with the following {str(len(genomeList))} taxa...')
speciesNames = set()
for species in genomeList:
	print(species.name)
	speciesNames.add(species.name)

# Write data into output database
with outPath.open(mode='w', newline='') as outfile:
	for protein in ref.db.values():
		for taxa in protein.taxa:
			if taxa in speciesNames:
				outfile.write(protein.getEntry())
				break

print('Database species optimizer completed successfully.')






















