#!/usr/bin/env python3
# coding: utf-8

import decimal
import clyngor
from clyngor import solve
from operator import itemgetter
from python_scripts.EMalgo import vraissemblance, expectation, maximization

ASP="""\

input(F):- matrixreadspecies(nf,F).

speciesread(S,I,W):- matrixreadspecies(nf,Filematrix); readsstrainsrank(I,J,W); species(J,S).
speciesread(S,I):- speciesread(S,I,_).


% Subsets names (species) and  elements of the set (reads)
subset(S) :- speciesread(S,_).
element(X) :- speciesread(_,X).

nb(subset(N)):- N={subset(X)}.
nb(element(N)):- N={element(X)}.
nb(weight(N)):- N=#sum{W,S,X:speciesread(S,X,W)}.


% Compute the threshold for the size of marginal subsets (under this threshold, subsets are considered marginal): 
% it's 1% of elements (reads)
onepercentelement((N+50)/100):- nb(element(N)).
% Compute the threshold for the weight of marginal subsets (under this threshold, subsets are considered marginal): 
% it's 2% of the total weight of elements (read ranks)
twopercentweight((N+25)/50):- nb(weight(N)).

marginal(S):- onepercentelement(N); subset(S); {speciesread(S,X)}N.
marginal(S):- twopercentweight(N); subset(S); #sum{W,X:speciesread(S,X,W)}N.
marginal_element(X):- marginal(S); speciesread(S,X).


% Choose the subsets, at least one
% The solution is made of non marginal selected strains
%---------------------------------------------------------
1{ choice(S):subset(S), not marginal(S) }.


%Statistics on the number/weight of elements covered by a chosen subset
wsupport(S,N):- choice(S); N=#sum{W,X:speciesread(S,X,W)}.


% A chosen subset (species) S is signed with an element (read) X 
%1) if it is the sole chosen subset covering this element
signedsubset1(S,X,W):- choice(S); speciesread(S,X,W); W>0; {speciesread(T,X):choice(T)}1.
%2) if there are exactly two subsets covering the element, but S has a better support than the other, which is not a signedsubset1 subset
%signedsubset2(S,X,W):- choice(S); speciesread(S,X,W); choice(S2); not signedsubset1(S2,_,_); speciesread(S2,X,W2); W>W2;  2{speciesread(T,X):choice(T)}2; not sol(S).
signedsubset2(S,X,W):- choice(S); speciesread(S,X,W); choice(S2); not signedsubset1(S2,_,_); speciesread(S2,X,W2); W>W2; W>1; 2{speciesread(T,X):choice(T)}2; not sol(S).

signedsubset(S,X,W):- signedsubset1(S,X,W).
signedsubset(S,X,W):- signedsubset2(S,X,W).
signedsubset(S,X):- signedsubset(S,X,_).


%Statistics on the number/weight of elements covered by an chosen signed subset
specwsupport(S,N):- N=#sum{W,X:signedsubset(S,X,W)}; N>0; choice(S). 


%An element (read) is specific if it is covered by at most one chosen subset
specific(X) :- element(X); { choice(S) : speciesread(S,X) } 1.

%Constraints
%-----------

% For each element, some subset must be chosen
:- element(X); not marginal_element(X); not choice(S) : speciesread(S,X).

% Each selection must be justified by some specific element
:- choice(S), not specific(X) : speciesread(S,X).


% Optimization
%-------------


%Minimize the number of chosen subsets that are not marginal
#minimize{1@10,S:choice(S)}.

%Maximize the total weighted support, then the total weighted specific support
#maximize {N@5,S:wsupport(S,N)}.
#maximize {N@3,S:specwsupport(S,N)}.


%Solution restricted to signed subsets
identified(S):- choice(S); specwsupport(S,N).


%Output
%------

#show parameters/3.
#show nb/1.
#show identified/1.
"""

def get_length(lengthFile: str) -> dict:
    '''Extract the length of each strains from the file'''
    
    lengths = {}
    with open(lengthFile,'r') as f:
        for line in f:
            line = line.strip().split('\t')
            genome = line[0].replace(".","_")
            length = line[1]
            lengths[genome]=length
    return lengths
            
def create_non_decorated_matrix(filename: str, outfilename: str):
    ''' Pierre: create the input_matrix removing headers (first line and first column)
    '''
    
    with open(filename, 'r') as fd, open(outfilename, 'w') as fd_out:
        fd.readline() # Pierre: do not consider first line
        for line in fd:
            line=line.strip().split()[1:]
            formated_line = " ".join(line)
            fd_out.write(f"{formated_line}\n")

            
        
def read_species_names(filename: str) -> str:
    '''Get the names of the strains.
    Return each atom species({number},"{name}")
    
    Modified by Pierre Peterlongo, input is now the kmindex output matrix'''
    
    names = []
    with open(filename) as fd:
        # read first line: 
        # test_coli       head_042_GCA_000027125_1_fasta  head_101-1_GCA_000168095_1_fasta        head_11368_GCA_000091005_1_fasta        head_12009_GCA_000010745_1_fasta        head_2H-327-20_GCA_902849165_1_fasta    head_4051-6_GCA_012532675_1_fasta  head_4608-58_GCA_000805835_1_fasta      head_53638_GCA_000167915_2_fasta        head_536_GCA_000013305_1_fasta  head_55989_GCA_000026245_1_fasta        head_921A_GCA_900499975_1_fasta    head_94389_GCA_001514625_1_fasta        head_APECO78_GCA_000332755_1_fasta      head_ATCC35469T_GCA_000026225_1_fasta   head_B1147_GCA_001660175_1_fasta        head_B253_GCA_000190495_1_fasta head_B7A_GCA_000725265_1_fasta     head_CFT073_GCA_000007445_1_fasta       head_CIP61_11_GCA_900236115_1_fasta     head_E1118_GCA_002109985_1_fasta        head_E24377A_GCA_000017745_1_fasta      head_ECOR01_GCA_002190105_1_fasta  head_ECOR24_GCA_002190595_1_fasta       head_ECOR37_GCA_002190835_1_fasta       head_ECOR42_GCA_002190935_1_fasta       head_ECOR49_GCA_002190975_1_fasta       head_ECOR60_GCA_002189835_1_fasta       head_ECOR63_GCA_002189905_1_fasta  head_ECOR64_GCA_002189945_1_fasta       head_ED1a_GCA_000026305_1_fasta head_FN-B26_GCA_902505495_1_fasta       head_H1-003-0072_GCA_902709615_1_fasta  head_H1-006-0003_GCA_902711215_1_fasta     head_H1-006-0025_GCA_902711405_1_fasta  head_H10407_GCA_000210475_1_fasta       head_H176_GCA_902842355_1_fasta head_H223_GCA_002110555_1_fasta head_H299_GCA_000176695_2_fasta head_HIPH08472_GCA_001514985_1_fasta       head_HS_GCA_000017765_1_fasta   head_NA114_GCA_000214765_3_fasta        head_RDEx444_GCA_003123505_1_fasta      head_S286_GCA_013363015_1_fasta head_S88_GCA_000026285_2_fasta  head_SE15_GCA_000010485_1_fasta    head_SMS-3-5_GCA_000019645_1_fasta      head_TW09231_GCA_000208465_2_fasta      head_UMN026_GCA_000026325_2_fasta
        line = fd.readline().strip().split()
        # remove first entry
        line = line[1:]
        for num, name in enumerate(line):
            name=name.rstrip("_fasta").rstrip("_fna").rstrip("_fq").rstrip("_fa").rstrip("_fastq")
            names.append(f'species({num},"{name}").')
    return('\n'.join(names)+'\n')

def read_matrix(filename: str,nbchoices: int,threshold: float) -> str:
    '''Get the strains for each read from the matrix. 
    Return each atom readsstrainsrank(read_number, strain_number, rank).''' 
    
    decimal.getcontext().prec = 2
    l=[]
    nbchoices = int(nbchoices)
    
    with open(filename) as fd:
        fd.readline() # Pierre: do not consider first line
        for i, line in enumerate(fd):
            sline=line.strip().split()[1:] # remove the read name (Pierre Peterlongo using kmindex output)
            vector = []
            for j,w in enumerate(sline):
                if (decimal.Decimal(w)*100) >= int(threshold):
                    vector.append((j+1,int(decimal.Decimal(w)*100))) 
            if len(vector) < int(nbchoices*1.5):
                # print(f"PIERRE (nbchoices,i,vector) {(nbchoices,i,vector)}")
                a = rank_list(nbchoices,i,vector)
                # print(f"PIERRE a is {a}")
                l.extend(a)
    return('\n'.join(l)+'\n')

def rank_list(n: int, ident: int, x) -> list[str]:
    ''' Gives the rank of at most n greatest elements greater than th in a list x with possible ties
    Returns triplets (id, index in x, rank in x)'''
    m = len(x)
    mn = min(n,m)
    s = sorted( x,  key=itemgetter(1), reverse=True)
    # print(f"PIERRE sorted s {s} of len {len(s)}")
    # print(f"PIERRE mn is {mn}")
    # ## DEBUG
    # s = [(1,100),(3,83), (5,83)]
    # print(f"PIERRE sorted s {s}")
    # mn=len(s)
    # ## ENDEBUG

    t = list(range(mn))
    for k in range(mn-1):
        t[k+1] = t[k] + (s[k+1][1] != s[k][1])
    
    # print(f"PIERRE strange t {t}")

    res = []
    for i,(j,_) in enumerate(s[:n]):
        res.append(f'readsstrainsrank({ident},{j},{t[mn-1]-t[i]}).')

    return res
    
    

def treat_strain(strain,dic):
    '''Take a dictionary and a set of strains and add 1 to the set of strains'''
    if frozenset(strain) in dic.keys() :
        dic[frozenset(strain)] += 1
    else:
        dic[frozenset(strain)] = 1

def quantification(file,strains,marginals):
    dic ={"Unknow":0,"Possibly marginal":0,"Possibly error":0}
    with open(file,'r') as f:
        f.readline()
        strain = set()
        marginal = set()
        for line in f:
            line = line.strip()
            if line[0] == '*':
                if line.split(' ')[1] == '0':
                    dic['Unknow'] += 1
                elif strain == set():
                    if marginal == set():
                        dic['Possibly error'] += 1
                    else:                        
                        dic['Possibly marginal'] += 1
                else:
                    treat_strain(strain,dic)
                strain = set()
                marginal = set()
            else:
                s = line.split(' ')[0].split('.fasta')[0]
                if s in strains:
                    strain.add(s)
                elif s in marginals:
                    marginal.add(s)
        if line.split(' ')[1] == '0':
            dic['Unknow'] += 1
        elif strain == set():
            if marginal == set():
                dic['Possibly error'] += 1
            else:                        
                dic['Possibly marginal'] += 1
        else:
            treat_strain(strain,dic)
    #print(dic)
    
def compute_abundances(out,strains,nb_strains,abondances,reads,length):
    '''Computation of the abondance, this part of the code still needs to be tested and improved'''
    
    print(f'Start : {abondances}')
   
    v = vraissemblance(strains,nb_strains,reads,abondances,length)
    
    for i in range(100):
        
        #Expectation part
        njs = expectation(strains,nb_strains,reads,abondances)

        #Maximisation part
        ajs = maximization(strains,nb_strains,reads,njs,length)
        
        modif = 0
        for i in range(nb_strains):
            modif = modif+(abondances[i]-ajs[i])
        abondances = ajs
        
        v = vraissemblance(strains,nb_strains,reads,abondances,length)
    
    #Write the results
    with open(out,'w') as o:
        for i,s in enumerate(strains):
            o.write(f'{s} : {round(abondances[i],2)}\n')

def main(args):
    
    clyngor.CLINGO_BIN_PATH = args.clingo_path
    
    #Get the strains names 
    names = read_species_names(args.matrix)
    print(f"Pierre names={names}")
    #Get the strains for every reads
    strainsreads = read_matrix(args.matrix,args.nbchoices,args.threshold)
    print(f"Pierre strainsreads={strainsreads}")
    
    #Get genomes length
    # here genome names still contain '.' symbols. 
    lengths = get_length(args.length)
    print(f"Pierre lenghts {lengths}")
    
    # Pierre: create a non decodated matrix in 
    outfilename="temp_matrix.tsv"
    create_non_decorated_matrix(args.matrix, outfilename)

    #Initialisation of the constant for the ASP script
    constant = f'#const nf=1. \n#const threshold={args.threshold}. \n#const nbchoices={args.nbchoices}.'
    matrix = f'matrixreadspecies(1,"{outfilename}").'
    
    #order the ASP script to run in toexec and solve it
    toexec = f"{constant}{matrix}{names}{strainsreads+ASP}"
    print(f"Pierre: toexec is {toexec}")
    answers = solve(inline=toexec,options='--opt-mode=optN')
    print(f"Pierre: SOLVED {answers}")
    for answer in answers: 
        print(f"Pierre answer {answer}")
    
    strains_in = set()
    i=0
    models = []
    for answer,optimization,optimality,answerNumber in answers.with_answer_number:    
        print(f'Answer {i+1}',optimization,optimality,answerNumber)
        if optimality:
            models.append(answer)
        i += 1
    
    #If no answer then stop the execution and ask for more reads
    if i == 0:
        print('No answer ! Maybe not enough reads or the stub (s) from the reads is not in the index.')
        exit(1)
    
    soluce = set()

    for model in models:
        for atom in model:
            if atom[0] == 'identified':
                strain = atom[1][0][1:-1]
                if '.bf' in strain:
                    strain = strain.split('.bf')[0]
                if '.fasta' in strain:
                    strain = strain.split('.fasta')[0]
                elif '.fna':
                    strain = strain.split('.fna')[0]
                strains_in.add(strain)
                
        soluce.add(frozenset(strains_in))

    strains_in = list(soluce)[-1]
    
    #In case of multiple optimum we take the set with the less strains
    for s in soluce:
        if len(s) <= len(strains_in):
            strains_in = set(s)
    #Here we have the strains present in the sample (strains_in)
    
    reads = []
    with open(args.file,'r') as f: #We read the results file from HowDeSBT
        f.readline()
        genomes = set()
        for line in f:
            if line[0] == '*':
                reads.append(genomes)
                genomes = set()
            else:
                line = line.strip().split(' ')
                if '.fna' in line[0] or '.fasta' in line[0]:
                    genome = ".".join(line[0].split('.')[:-1])
                else:
                    genome= line[0]
                if '.fna' in line[0] or '.fasta' in line[0]:
                    if genome[-1] == '_':
                        genome = genome[:-1]
                if genome in strains_in:
                    genomes.add(genome)
    reads.append(genomes)
    
    strains = sorted(list(strains_in))
    nb_strains = int(len(strains))
    
    #Start the abondance computation with all strains at the same level
    abondances = []
    for i in range(nb_strains):
        abondances.append(float(1/nb_strains))
        
    #Get the length of each strains from the sample
    strain_lengths = []
    for s in strains:
        try :
            strain_lengths.append(float(lengths[s]))
        except KeyError:
            s = '-'.join(s.split('_'))
            strain_lengths.append(float(lengths[s]))
    
    #Compute the abondance for each strain and write the output file
    compute_abundances(args.output,strains,nb_strains,abondances,reads,strain_lengths)
        

if __name__ == "__main__":
    main()


