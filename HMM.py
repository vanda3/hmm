from itertools import product
import csv
from graphviz import Digraph
import math
import numpy as np
import sys

sys.setrecursionlimit(100000)

debug=False
v=np.empty(2)
f=np.empty(2)
trace=np.empty(2)

class state:
    def __init__(self, name, bases, idx):
        self.idx=idx
        self.name=name
        self.keys=bases.keys()
        self.bases=bases
        self.total_e=0
        self.total_t=0
        self.next=[]
        self.transition={}
        self.prob_t={}
        self.prob_e={}
        self.prob_start={}
    def e(self,name):
        self.bases[name]+=1
        self.total_e+=1
    def t(self, name):
        self.transition[name]+=1
        self.total_t+=1
    def edge(self,nxt):
        self.transition[nxt.name]=0
        self.next.append(nxt)


################################################################################
# AUX
################################################################################
def calculate(states):
    global debug
    # emissions
    for s in states:
        for k in states[s].bases:
            if(states[s].total_e>0 and states[s].bases[k]>0):
                states[s].prob_e[k]=states[s].bases[k]/states[s].total_e
            else:
                states[s].prob_e[k]=0.001
            if(debug):
                print("state:", s, " e:",k, " prob=",states[s].prob_e[k])
    # transitions
    for s in states:
        for k in states[s].transition:
            if(states[s].total_t>0 and states[s].transition[k]>0):
                states[s].prob_t[k]=states[s].transition[k]/states[s].total_t
            else:
                states[s].prob_t[k]=0.001
            if(debug):
                print(s, " -> ",k, " prob=",states[s].prob_t[k])
    states["sNG"].prob_e["sNG"]=1
    states["fNG"].prob_e["fNG"]=1
    states["sG"].prob_e["sG"]=1
    states["fG"].prob_e["fG"]=1
    states["start"].prob_e["start"]=1


################################################################################
# AUX
################################################################################
def countEvents(seq,start,stop,direction):
    enterNG=False
    for i in range(0,len(start)):
        if(debug):
            print("////// i:",i)
        if i==0 and start[0]>0:
            states["start"].t("sNG")
            b=0
            e=start[0]
            if e>b:
                if(debug):
                    print("b:",b,", e:",e)
                nongene(seq[b:e],states)
                states["fNG"].t("sG")
                enterNG=True
                if(debug):
                    print("start:",start[i],", stop:",stop[i], " direction:",direction[i])
        if(i==0 and start[0]<=0):
            states["start"].t("sG")
        if(direction[i]=="+"):
            gene(seq[start[i]:stop[i]+1],states)
            enterNG=False
            if(debug):
                print("GENE: ",seq[start[i]:stop[i]+1])
        else:
            gene(reverse_dna(seq[start[i]:stop[i]+1]),states)
            enterNG=False
            if(debug):
                print("GENE INVERSO: ",reverse_dna(seq[start[i]:stop[i]+1]))
        states["fG"].t("sNG")
        if i!=len(start)-1:
            b=stop[i]+1
            e=start[i+1]
            if (e-1)>b:
                if(debug):
                    print("b:",b,", e:",e)
                nongene(seq[b:e],states)
                states["fNG"].t("sG")
                enterNG=True
        if i==len(start)-1:
            if(stop[len(start)-1]<len(seq)):
                b=stop[len(start)-1]+1
                if(b<len(seq)-1):
                    if(debug):
                        print("b:",b,", e:",len(seq))
                    nongene(seq[b:len(seq)],states)
                    states["fNG"].t("sG")
                    enterNG=True
        if(enterNG==False):
            states["fG"].t("sG")


################################################################################
# PRINT
################################################################################
def printHMM(g, states):
    for s in states:
        for k in states[s].next:
            if(states[s].total_t>0):
                g.edge(states[s].name, k.name, label=str(round(states[s].prob_t[k.name],2)))


################################################################################
# PARSE
################################################################################
def parse(name):
    name=str(name)+".txt"
    file = open(name,mode='r+')
    sequence=""
    for line in file:
        line=line.strip()
        sequence+=line
    return sequence


################################################################################
# PARSECSV
################################################################################
def parseCSV(name):
    start=[]
    stop=[]
    direction=[]
    i=0
    with open(name+'.csv', 'r+') as f:
        reader = csv.reader(f)
        next(reader, None)
        for row in reader:
            row=row[0].split(",")
            start.append(int(row[0])-1)
            stop.append(int(row[1])-1)
            direction.append(str(row[2]).strip('"'))
            i+=1
    return start, stop, direction

################################################################################
# AUX PARSE
################################################################################
def gene(seq,states):
    st=["e1","e2","e3","e4","e5"]
    curr=states["sG"]
    if(seq[0:2]=="GC"):
        return
    for i in range(0,len(seq)-3):
        char=seq[i]
        for s in curr.next:
            if((char in s.keys) and (s.name not in st)):
                states[curr.name].t(s.name)
                states[s.name].e(char)
                #print(curr.name, " -",char,"-> ",s.name)
                curr=states[s.name]
    char=seq[i+1]
    states[curr.name].t("e1")
    #print(curr.name, " -",char,"-> e1")
    curr=states["e1"]
    states["e1"].e(char)
    for i in range(i+2,len(seq)):
        char=seq[i]
        for s in curr.next:
            if(char in s.keys):
                states[curr.name].t(s.name)
                states[s.name].e(char)
                #print(curr.name, " -",char,"-> ",s.name)
                curr=states[s.name]
    if(s.name=="e4"):
        states["e4"].t("fG")
    if(s.name=="e5"):
        states["e5"].t("fG")


################################################################################
# AUX PARSE
################################################################################
def nongene(seq,states):
    global debug

    if(debug):
        print("NON GENE: ",seq)
    if(seq[0]=="A"):
        states["sNG"].t("A")
    elif(seq[0]=="C"):
        states["sNG"].t("C")
    elif(seq[0]=="G"):
        states["sNG"].t("G")
    else:
        states["sNG"].t("T")
    states[seq[0]].e(seq[0])
    char=""
    char2=""
    for i in range(0,len(seq)-1):
        char=seq[i]
        char2=seq[i+1]
        states[char].t(char2)
        states[char2].e(char2)
    states[char2].t("fNG")


################################################################################
# AUX PARSE
################################################################################    
def reverse_dna(sub_string):
	associate = {"A":"T","C":"G","T":"A","G":"C"}
	reverse_string = sub_string[::-1]
	dna =""
	for c in reverse_string:
		dna += associate[c]
	return dna


################################################################################
# AUX MAIN
################################################################################
def build():
    states={}
    ls={}
    b = open("build.txt", 'r')
    i=0
    idx=0
    for line in b:
        if i==0:
            n_states=int(line.strip())
        elif i in range(1,n_states+1):
            l=line.split()
            bases={}
            name=l[0].strip()
            for j in range(1,len(l)):
                bases[l[j].strip()]=0
            states[name]=state(name,bases,idx)
            ls[idx]=states[name]
            idx+=1
        elif i ==(n_states+1):
            n_edges=int(line.strip())
        else:
            e=line.split()
            states[e[0].strip()].edge(states[e[1].strip()])
        i+=1
    b.close()
    return states,ls


################################################################################
# VITERBI
################################################################################
def viterbi(s,seq,training):
    global debug
    global v
    global trace

    lines=len(s)
    cols=len(seq)
    print("//// LINES: ",lines," COLS:", cols," ////")
    v=np.full((lines, cols),0.0)
    trace=np.full((lines, cols),-1)
    for k in s[0].next: #sNG e sG
        for j in k.next:
            if(seq[0] in j.bases):
                if(debug):
                    print("/////////// start -> ",k.name," -",seq[0],"-> ",j.name)
                v[j.idx][0]=math.log(s[0].prob_t[k.name],2)+math.log(k.prob_t[j.name],2)+math.log(j.prob_e[seq[0]],2)
                if(not training):
                    trace[j.idx][0]=0
                recursiveV(s,seq,j.idx,1,lines,cols,training)
    if(not training):
        print("Probability Matrix:")  
        printMatrix(v,lines,cols)
        print("Trace Matrix:")
        printMatrix(trace,lines,cols)
        traceback(s,lines,cols)
    if(training):
        train(s,seq,"v",lines,cols)


################################################################################
# AUX VITERBI
################################################################################
def recursiveV(s,seq,prev_line,col,lines,cols,training):
    global v
    global debug
    global trace
    
    b=["A","C","G","T"]
    if(col==cols):
        return
    for i in range(1,lines):
        if(s[i] in s[prev_line].next and s[i]=="fNG"):
            j=s[i].next[0]
            for k in j.next:
                if(seq[col] in k.bases):
                    if(debug):
                        print("/////////// ",s[prev_line].name," -> fNG -> ",s[i].name," -",seq[col],"-> ",k.name)
                    v[k.idx][col]=v[prev_line][col-1]+math.log(s[prev_line].prob_t["fNG"],2)+math.log(j.prob_t[k.name],2)+math.log(k.prob_e[seq[col]],2)
                    trace[k.idx][col]=prev_line
                    recursiveV(s,seq,k.idx,col+1,lines,cols,training)
        elif(s[i] in s[prev_line].next and s[i]=="fG"):
            for j in s[i].next:
                if(j.name=="sNG"):
                    for k in j.next:
                        if(seq[col] in k.bases):
                            if(debug):
                                print("/////////// ",s[prev_line].name," -> fG -> sNG -",seq[col],"-> ",k.name)
                            v[k.idx][col]=v[prev_line][col-1]+math.log(s[prev_line].prob_t[s[i].name],2)+math.log(s[i].prob_t[j.name],2)+math.log(j.prob_t[k.name],2)
                            trace[k.idx][col]=prev_line
                            recursiveV(s,seq,k.idx,col+1,lines,cols,training)
                else:
                    for k in j.next:
                        if(seq[col] in k.bases):
                            if(debug):
                                print("/////////// ",s[prev_line].name," -> fG -> sG -",seq[col],"-> ",k.name)
                            v[k.idx][col]=v[prev_line][col-1]+math.log(s[prev_line].prob_t[s[i].name],2)+math.log(s[i].prob_t[j.name],2)+math.log(j.prob_t[k.name],2)+math.log(k.prob_e[seq[col]],2)
                            trace[k.idx][col]=prev_line
                            recursiveV(s,seq,k.idx,col+1,lines,cols,training)          
        elif(s[i] in s[prev_line].next and s[i] in b and seq[col] in s[i].bases):
            if(debug):
                print("/////////// ",s[prev_line].name," -",seq[col],"-> ",s[i].name)
            v[s[i].idx][col]=v[prev_line][col-1]+math.log(s[prev_line].prob_t[s[i].name],2)
            trace[s[i].idx][col]=prev_line
            recursiveV(s,seq,s[i].idx,col+1,lines,cols,training)
        elif(s[i] in s[prev_line].next and seq[col] in s[i].bases):
            if(debug):
                print("/////////// ",s[prev_line].name," -",seq[col],"-> ",s[i].name)
            v[s[i].idx][col]=v[prev_line][col-1]+math.log(s[prev_line].prob_t[s[i].name],2)+math.log(s[i].prob_e[seq[col]],2)
            trace[s[i].idx][col]=prev_line
            recursiveV(s,seq,s[i].idx,col+1,lines,cols,training)
        else:
            continue


################################################################################
# AUX VITERBI
################################################################################
def traceback(s,lines,cols):
    col=cols-1
    path=[]
    maxi=-1000000
    pos=0
    for i in range(0,lines):
        if v[i][col]>maxi and v[i][col]!=0.0:
            maxi=v[i][col]
            pos=i
    path.append(s[pos].name)
    prev_state=trace[pos][col]
    path.append(s[prev_state].name)
    col-=1
    while(col>0):
        prev_state=trace[prev_state][col]
        path.append(s[prev_state].name)
        col-=1
    path.append("start")
    k=len(path)
    print("PATH: ",end='')
    while(k>0):
        k-=1
        print(path[k]," ",end='')
    print()


################################################################################
# AUX VITERBI
################################################################################
def printMatrix(matrix,lines,cols):
    for i in range(0,lines):
        for j in range(0,cols):
            print(matrix[i][j]," ",end='')
        print()  


################################################################################
# FORWARD
################################################################################
def forward(states,s,seq,training):
    global debug
    global f
    global trace

    lines=len(s)
    cols=len(seq)
    if(debug):
        print("//// LINES: ",lines," COLS:", cols," ////")
    f=np.full((lines, cols),0.0)
    trace=np.full((lines, cols),-1)

    for k in s[0].next: #sNG e sG
        for j in k.next:
            if(seq[0] in j.bases):
                if(debug):
                    print("/////////// start -> ",k.name," -",seq[0],"-> ",j.name)
                f[j.idx][0]=math.log(s[0].prob_t[k.name],2)+math.log(k.prob_t[j.name],2)+math.log(j.prob_e[seq[0]],2)
                trace[j.idx][0]=0

    for j in range(1,cols):
        for i in range(1,lines):
            result=recursiveF(states,s,seq,i,j,lines)
            if(result!=0):
                f[i][j]=result

    if(not training):
        print("Probability Matrix:")  
        printMatrix(f,lines,cols)

    if(training):
        train(s,seq,"f",lines,cols)



################################################################################
# AUX FORWARD
################################################################################
def recursiveF(states,s,seq,i,col,lines):
    global f
    global debug
    global trace

    soma=0
    b=["A","C","G","T"]
    nulls=["sG","fNG","fG","sNG"]
    for prev_line in range(1,lines):
        if((s[i].name not in nulls) and (s[prev_line].name not in nulls) and s[i]):
            if((s[i].name == "s1" or s[i].name=="s2") and (seq[col] in s[i].bases)):
                # e4/5 -> s1/2
                if(s[prev_line].name=="e4" or s[prev_line].name=="e5"):
                    if(debug):
                        print("/////////// ",s[prev_line].name," fG -> sG -",seq[col],"-> ",s[i].name, " i:",i," col:",col)
                    soma+=f[prev_line][col-1]+math.log(s[prev_line].prob_t["fG"],2)+math.log(states["fG"].prob_t["sG"],2)+math.log(states["sG"].prob_t[s[i].name],2)+math.log(s[i].prob_e[seq[col]],2)
                    trace[i][col]=prev_line
                # A/C/T/G -> s1/2
                elif(s[prev_line].name in b):
                    if(debug):
                        print("/////////// ",s[prev_line].name," fNG -> sG -",seq[col],"-> ",s[i].name)
                    soma+=f[prev_line][col-1]+math.log(s[prev_line].prob_t["fNG"],2)+math.log(states["fNG"].prob_t["sG"],2)+math.log(states["sG"].prob_t[s[i].name],2)
                    trace[s[i].idx][col]=prev_line

            elif((s[i].name in b) and (seq[col] in s[i].bases)):
                # A/C/T/G -> A/C/T/G
                if(s[prev_line] in b):
                    if(debug):
                        print("/////////// ",s[prev_line].name," -",seq[col],"-> ",s[i].name)
                    soma+=f[prev_line][col-1]+math.log(s[prev_line].prob_t[s[i].name],2)
                    trace[s[i].idx][col]=prev_line
   
               ## e4/5 -> A/C/T/G
                if(s[prev_line]=="e4" or s[prev_line]=="e5"):
                    if(debug):
                        print("/////////// ",s[prev_line].name," fG -> sNG -",seq[col],"-> ",s[i].name)
                    soma+=f[prev_line][col-1]+math.log(s[prev_line].prob_t["fG"],2)+math.log(states["sNG"].prob_t["sNG"],2)+math.log(states["sG"].prob_t[s[i].name],2)
                    trace[s[i].idx][col]=prev_line

            elif(s[i] in s[prev_line].next and (seq[col] in s[i].bases)):
                if(debug):
                    print("/////////// ",s[prev_line].name," -",seq[col],"-> ",s[i].name)
                soma+=f[prev_line][col-1]+math.log(s[prev_line].prob_t[s[i].name],2)+math.log(s[i].prob_e[seq[col]],2)
                trace[s[i].idx][col]=prev_line
            
            else:
                continue
    return soma


################################################################################
# AUX ALGS
################################################################################
def train(s,seq,alg,lines,cols):
    global v
    global f
    global trace

    A=np.full((lines, cols),0.0)
    E=np.full((lines, cols),0.0)

    pseq=0

    if(alg=="v"):
        m=v
    else:
        m=f

    for i in range(1,lines):
        pseq+=f[i][cols-1]

    # update A
    for k in range(0,lines-1):
        for l in range(0,lines-1):
            esq=1/pseq
            dire=0
            for i in range(0,cols):
                if(trace[k][l]!=-1 and seq[i] in s[l].bases):
                    dire+=m[k+1][i]+math.log(s[k].prob_t[s[l].name],2)+math.log(s[l].prob_e[seq[i]],2)
                if(trace[k][l]!=-1 and seq[i]==s[l].name):
                    dire+=m[k+1][i]+math.log(s[k].prob_t[s[l].name],2)
            A[k][l]+=esq+dire

    # update E
    for k in range(0,lines-1):
        for l in s[k].bases: 
            esq=1/pseq
            dire=0
            for i in range(0,cols):
                if(trace[k][l]!=-1 and seq[i-1] == l):
                    dire+=m[k+1][i]
            E[k][l]+=esq+dire


################################################################################
# MAIN
################################################################################
if __name__ == '__main__':
    i=0
    order=1 
    stop=[]
    start=[]
    direction=[]
    states,l_states=build()
    name="ecolis"

    seq=parse(name)
    start,stop,direction=parseCSV(name)
    countEvents(seq,start,stop,direction)
    calculate(states)
    
    g = Digraph('G', filename='HMM')
    printHMM(g, states)
    g.view()

    print("Insert sequence:")
    seq=str(input())

    print("Choose training algorithm:")
    print("1) Viterbi")
    print("2) Forward")

    opt=int(input())
    if(opt==1):
        viterbi(l_states,seq,False)
    else:
        forward(states,l_states,seq,False)