from matplotlib import pyplot as plt
import os, pickle, math, sys, random
import numpy as np
from scipy.stats import norm
import pandas as pd

class Stats():
    def __init__(self,path):
        if path != '' and path[-1] != '/':
            path += '/'
        file = open(path+'component_stats.txt')
        f = file.read()
        file.close()
        f = f.strip().splitlines()
        self.stats = [line.split('\t')[0] for line in f]
        self.stat2score = {s:np.nan for s in self.stats}
    def set_stat(self,stat,value):
        self.stat2score[stat] = value


class Mix1D():
    def __init__(self):
        self.means_ = []
        self.weights_ = []
        self.covars_ = []

class ODE_HMM:
    '''
    emissions use ODE terms -- imputation is done based on "ode_compensation"
    '''

    def __init__(self, path2trained):
        self._initialize_AODE(path2trained)
        self.forcesweep=True
        self.hmmstate2trainingstate = {'Neutral':'neutral','LinkedLeft':'linked','LinkedRight':'linked','Sweep':'sweep'}
        self.hmmstate2num = {'StartEnd':0,'Neutral':1,'LinkedLeft':2,'LinkedRight':3,'Sweep':4}
        self.num2hmmstate = {y:x for x,y in self.hmmstate2num.items()}
        self.hmmstates = self.hmmstate2num.keys()
        self.set_transitions()


    def set_transitions(self):
        self.transitions = [[0 for j in range(len(self.hmmstates))] for i in range(len(self.hmmstates))]

        #Start/End goes to Neutral (1) or LinkedLeft (2)
        self.transitions[0][1] = 0.99999
        self.transitions[0][2] = 0.00001

        #Neutral goes to Neutral (1) or LinkedLeft (2) or StartEnd (0)
        # self.transitions[1][1] = 0.999998
        # self.transitions[1][2] = 0.000001
        # self.transitions[1][0] = 0.000001

        self.transitions[1][1] = 0.999998999
        self.transitions[1][2] = 0.000000001
        self.transitions[1][0] = 0.000001       

        if self.forcesweep:
        #LinkedLeft goes to LinkedLeft (2) or Sweep (4)
            self.transitions[2][2] = 0.999
            self.transitions[2][4] = 0.001

        else: #allow back to neutral
            self.transitions[2][2] = 0.998
            self.transitions[2][4] = 0.001
            self.transitions[2][1] = 0.002

        #Sweep goes to LinkedRight (3)
        self.transitions[4][3] = 1

        #LinkedRight goes to LinkedRight (3) or Neutral (1)
        self.transitions[3][3] = 0.999
        self.transitions[3][1] = 0.001  

    def logsumexp(self,sumvec): #sumvec is log of entries in the sum (x_i's in log-sum-exp)
        #e.g. logV[i-1][j']+log(t[j'][j])
        a = max(sumvec)
        exponents = [x-a for x in sumvec]
        newsum = sum([math.exp(x) for x in exponents])
        newlogsum = math.log(newsum)
        return a+newlogsum


    def viterbi(self,Svec,positions,outpath,outname,plotpaths=False,classify=True): #Svec[i] is an instance of Stats() for position i
        #V = [[0 for x in self.hmmstates] for i in range(len(Svec)+1)] #V[i][j] is for position i and state j
        logV = [['n/a' for x in self.hmmstates] for i in range(len(Svec)+1)]
        P = [['n/a' for x in self.hmmstates] for i in range(len(Svec)+1)] #P[i][j] is the state pointer from position i and state j to state at i-1
        #initialize (0 (or n/a) unless set):
        
        logV[0][1] = math.log(self.transitions[self.hmmstate2num['StartEnd']][self.hmmstate2num['Neutral']]) + math.log(self.emission('Neutral',Svec[0])[0])
        logV[0][2] = math.log(self.transitions[self.hmmstate2num['StartEnd']][self.hmmstate2num['LinkedLeft']]) + math.log(self.emission('LinkedLeft',Svec[0])[0])

        P[0][1] = 0 #pointer at position 0 and state Neutral (1) is to StartEnd (0)
        P[0][2] = 0 #pointer at position 0 and state LinkedLeft (2) is to StartEnd (0)
        
        print('Viterbi recursion...')
        #recursion
        
        for i in range(1,len(Svec)):
            #print positions[i]
            for j in range(1,len(self.hmmstates)):
                maxpathprobs = -1*np.inf
                for jprime in range(len(self.hmmstates)):
                    if logV[i-1][jprime] != 'n/a' and self.transitions[jprime][j] != 0:
                        logpp = logV[i-1][jprime] + math.log(self.transitions[jprime][j])
                        if logpp > maxpathprobs:
                            maxpathprobs = logpp
                            argmax = jprime
                logV[i][j] = maxpathprobs + math.log(self.emission(self.num2hmmstate[j],Svec[i])[0])
                P[i][j] = argmax

        #end state
        maxpathprobs = -1*np.inf
        for jprime in range(1,len(self.hmmstates)):
            if logV[len(Svec)-1][jprime] != 'n/a' and self.transitions[jprime][0] != 0:
                logpp = logV[len(Svec)-1][jprime] + math.log(self.transitions[jprime][0])
                if logpp > maxpathprobs:
                    maxpathprobs = logpp
                    argmax = jprime
        logV[len(Svec)][0] = maxpathprobs
        P[len(Svec)][0] = argmax 

        #viterbi path: follow pointers in P:
        state_path = [0 for i in range(len(Svec))]
        state_path[len(Svec)-1] = P[len(Svec)][0]
        for i in range(len(Svec)-2,-1,-1):
            state_path[i] = P[i+1][state_path[i+1]]

        if classify:
            out = open(opath.join(outpath, outname+'_viterbi_classified.txt'))
            out.write('pos\tstate\n')
            for i in range(len(Svec)):
                position = positions[i]
                out.write('\t'.join(str(position), state_path[i])+'\n')
            out.close()
    

        if plotpaths:
            plt.plot(positions,state_path,'o-',mfc='Salmon',mec='none')
            plt.savefig(outpath+outname+'_viterbi_path.pdf')
            plt.close()

    def backward(self,Svec):
        logV = [['n/a' for x in self.hmmstates] for i in range(len(Svec))]
        for state in self.hmmstates:
            if state != 'StartEnd' and self.transitions[self.hmmstate2num[state]][self.hmmstate2num['StartEnd']] != 0:
                logV[len(Svec)-1][self.hmmstate2num[state]] = math.log(self.transitions[self.hmmstate2num[state]][self.hmmstate2num['StartEnd']])
        print('Backward recursion...')
        for i in range(len(Svec)-2,-1,-1):
            for j in range(1,len(self.hmmstates)):
                sumvec = []
                for jprime in range(len(self.hmmstates)):
                    if logV[i+1][jprime] != 'n/a' and self.transitions[j][jprime] != 0:
                        sumvec.append(logV[i+1][jprime] + math.log(self.transitions[j][jprime]) + math.log(self.emission(self.num2hmmstate[jprime],Svec[i+1])[0]))
                if len(sumvec) > 0:
                    Sum = self.logsumexp(sumvec)
                    logV[i][j] = Sum
        Beginsumvec = []
        for jprime in range(len(self.hmmstates)):
            if logV[0][jprime] != 'n/a' and self.transitions[0][jprime] != 0:
                #print jprime
                Beginsumvec.append(logV[0][jprime]+math.log(self.transitions[0][jprime])+math.log(self.emission(self.num2hmmstate[jprime],Svec[0])[0]))
        BeginSum = self.logsumexp(Beginsumvec)


        return logV


    def forward(self,Svec):
        logV = [['n/a' for x in self.hmmstates] for i in range(len(Svec)+1)]
        logV[0][1] = math.log(self.transitions[self.hmmstate2num['StartEnd']][self.hmmstate2num['Neutral']]) + math.log(self.emission('Neutral',Svec[0])[0])
        logV[0][2] = math.log(self.transitions[self.hmmstate2num['StartEnd']][self.hmmstate2num['LinkedLeft']]) + math.log(self.emission('LinkedLeft',Svec[0])[0])

        print('Forward recursion...')   
        for i in range(1,len(Svec)):
            for j in range(1,len(self.hmmstates)):
                #Use: log(sum(exp(xn))) = a+log(sum(exp(xn-a))) where a=max_n(xn)               
                sumvec = []
                for jprime in range(len(self.hmmstates)):
                    if logV[i-1][jprime] != 'n/a' and self.transitions[jprime][j] != 0:
                        sumvec.append(logV[i-1][jprime]+math.log(self.transitions[jprime][j]))
                if len(sumvec) > 0:
                    Sum = self.logsumexp(sumvec)
                    logV[i][j] = Sum + math.log(self.emission(self.num2hmmstate[j],Svec[i])[0])
        Endsumvec = []
        for jprime in range(1,len(self.hmmstates)):
            if logV[len(Svec)-1][jprime] != 'n/a' and self.transitions[jprime][0] != 0:
                Endsumvec.append(logV[len(Svec)-1][jprime] + math.log(self.transitions[jprime][0]))
        EndSum = self.logsumexp(Endsumvec)
        logV[len(Svec)][0] = EndSum
        return logV                 

    def emission(self,scenario,S):

        scenario = self.hmmstate2trainingstate[scenario]
        Likelihoods = []
        for stat in self.statlist:
            L = self.ode_likelihood(stat,S,scenario)
            if L != 'n/a':
                Likelihoods.append(L)
        Likelihood = float(sum(Likelihoods))/len(Likelihoods) #average all ODE likelihoods -- note that missing statistics don't cause problems for this average.
        if Likelihood < 1e-200:
            Likelihood = 1e-200
        return Likelihood, Likelihoods

    def stochastic_backtrace(self,Svec,logV):
        #logV = self.forward(Svec)
        state_path = [0 for i in range(len(Svec))]
        #last position: logV[len(Svec)][0]
        unnormalized_logprobs = ['n/a' for i in range(len(self.hmmstates))]
        
        for jprime in range(len(self.hmmstates)):
            if logV[len(Svec)-1][jprime] != 'n/a' and self.transitions[jprime][0] != 0:
                unnormalized_logprobs[jprime] = logV[len(Svec)-1][jprime] + math.log(self.transitions[jprime][0])
        #print unnormalized_logprobs
        normalization = self.logsumexp([x for x in unnormalized_logprobs if x != 'n/a'])
        normalized_logprobs = ['n/a' for i in range(len(self.hmmstates))]
        for jprime in range(len(self.hmmstates)):
            if unnormalized_logprobs[jprime] != 'n/a' and self.transitions[jprime][0] != 0:
                normalized_logprobs[jprime] = logV[len(Svec)-1][jprime] + math.log(self.transitions[jprime][0])-normalization
        Probs = [0 for i in range(len(self.hmmstates))]
        for jprime in range(len(self.hmmstates)):
            if normalized_logprobs[jprime] != 'n/a':
                Probs[jprime] = math.exp(normalized_logprobs[jprime])
        #print Probs

        #Probs = [0 for i in range(len(self.hmmstates))]
        #for jprime in range(1,len(self.hmmstates)):
        #   Probs[jprime] = math.exp(logV[len(Svec)-1][jprime])*self.transitions[jprime][0]
        #S = sum(Probs)
        #Probs = [float(x)/S for x in Probs]
        state = np.random.choice(len(self.hmmstates),p=Probs)
        #print state
        state_path[len(Svec)-1] = state

        #print len(Svec)-2
        for i in range(len(Svec)-2,-1,-1):
            #print i
            unnormalized_logprobs = ['n/a' for thing in range(len(self.hmmstates))]
            for jprime in range(len(self.hmmstates)):
                if logV[i][jprime] != 'n/a' and self.transitions[jprime][state_path[i+1]] != 0:
                    unnormalized_logprobs[jprime] = logV[i][jprime] + math.log(self.transitions[jprime][state_path[i+1]])
            #if i == len(Svec)/2:
                #print unnormalized_logprobs
            normalization = self.logsumexp([x for x in unnormalized_logprobs if x != 'n/a'])
            normalized_logprobs = ['n/a' for thing in range(len(self.hmmstates))]
            for jprime in range(len(self.hmmstates)):
                if unnormalized_logprobs[jprime] != 'n/a' and self.transitions[jprime][state_path[i+1]]:
                    normalized_logprobs[jprime] = logV[i][jprime] + math.log(self.transitions[jprime][state_path[i+1]])-normalization

            
            Probs = [0 for thing in range(len(self.hmmstates))]
            for jprime in range(len(self.hmmstates)):
                if normalized_logprobs[jprime] != 'n/a':
                    Probs[jprime] = math.exp(normalized_logprobs[jprime])
            #print Probs

            # for jprime in range(1,len(self.hmmstates)):
            #   if logV[i][jprime] != 'n/a':
            #       Probs[jprime] = math.exp(logV[i][jprime])*self.transitions[jprime][state_path[i+1]]
            # S = sum(Probs)
            # Probs = [float(x)/S for x in Probs]
            state_path[i] = np.random.choice(len(self.hmmstates),p=Probs)
        #print state_path
        return state_path


    def many_backtraces(self,Svec,N,positions,outpath,outname,plotpaths=False,plotdensity=False,classify=True,plot_classify=True): ##this is the function to call to classify using the HMM
        '''
        Svec: list of statistic values across genomic region
        N: number of backtraces
        positions: phys positions to match Svec
        outpath: where to put output
        outname: how to name output (e.g. parameter_model_name)
        plotpaths (bool): plot all backtrace paths
        plotdensity (bool): plot frequency of sweep locations across all backtrace paths
        classify (bool): return classifications for each position based on backtrace paths
        plot_classify (bool): make a bar plot of classification proportions at each position (only use with classify=True)
        '''

        logV = self.forward(Svec)
        State_path_counts = [[0,0,0,0,0] for i in range(len(Svec))]
        #State_path_probabilities = [[0,0,0,0,0] for i in range(len(Svec))]
        totalsweeps = 0
        for i in range(N):
            statepath = self.stochastic_backtrace(Svec,logV)
            wentthroughsweep = False
            for j in range(len(statepath)):
                State_path_counts[j][statepath[j]] += 1
                if statepath[j] == 4:
                    wentthroughsweep = True
            if wentthroughsweep:
                totalsweeps += 1
            if plotpaths:
                plt.plot(positions,statepath,'o-')
        sweeppercent = float(totalsweeps)/N
        
        if plotpaths:
            plt.savefig(outpath+outname+'_backtrace_paths.pdf')
            plt.clf()


        # at each position, count fraction assigned to neutral/linked/sweep
        if classify:
            out = open(opath.join(outpath, outname+'_hmm_classified.txt'))
            out.write('pos\tP(neutral)\tP(linked)\tP(sweep)\n')
            n_prop_vec = []
            l_prop_vec = []
            s_prop_vec = []
            for i in range(len(Svec)):
                position = positions[i]
                n_count = State_path_counts[i][1]
                l_count = State_path_counts[i][2]+State_path_counts[i][3]
                s_count = State_path_counts[i][4]
                total = n_count+l_count+s_count
                n_prop = float(n_count)/total
                l_prop = float(l_count)/total
                s_prop = float(s_count)/total
                n_prop_vec.append(n_prop)
                l_prop_vec.append(l_prop)
                s_prop_vec.append(s_prop)
                out.write('\t'.join(str(position), n_prop, l_prop, s_prop)+'\n')
            out.close()
            if plot_classify:
                p1 = plt.bar(range(len(n_prop_vec)), n_prop_vec, width=0.8, color='blue')
                p2 = plt.bar(range(len(n_prop_vec)), l_prop_vec, width=0.8, bottom=n_prop_vec, color='green')
                p3 = plt.bar(range(len(n_prop_vec)), s_prop_vec, width=0.8, bottom=[n_prop_vec[i]+l_prop_vec[i] for i in range(len(n_prop_vec))], color='red')
                plt.xlabel('genomic position index')
                plt.ylabel('proportion of backtrace paths')
                plt.legend([p1[0], p2[0], p3[0]], ['Neutral', 'Linked', 'Sweep'])
                plt.savefig(opath.join(outpath, outname+'_hmm_classified.pdf'))
                plt.clf()



        # #concentrate on sweep state
        # V = [State_path_counts[i][4] for i in range(len(Svec))]
        # total = sum(V)
        # numpaths = float(total)/N
        # out = open(outpath+plotname+'_backtrace_density.txt','w')
        # for i in range(len(Svec)):
        #   if V[i] > 0:
        #       out.write(str(positions[i])+'\t'+str(V[i])+'\n')
        # out.close()
        if plotdensity:
            plt.plot(positions,V,color='maroon')
            plt.figtext(0.7,0.7,str(int(sweeppercent*100))+"% of paths\nhave a sweep")
            plt.xlabel('Genomic Position')
            plt.ylabel('Number of Paths')
            plt.savefig(outpath+outname+'_backtrace_density.pdf')
            plt.clf()

        return sweeppercent

    # def backtrace_density(self,datafile,outpath,plotname,N=10000): ##this is what to call to 
    #   Svec,positions = self.read_Svec(datafile,impute=self.impute)
    #   if len(Svec) > 0:
    #       sweeppercent = self.many_backtraces(Svec,N,positions,datafile,outpath,plotname,plotpaths=False,plotdensity=True)
    #   out = open(outpath+'sweep_percentages.txt','a')
    #   out.write(str(sweeppercent)+'\n')
    #   out.close()

####### read in allstats file and impute ###################

    def read_Svec(self,datafile):
        
        #S = Stats(self.path2trained)
        df = pd.read_csv(datafile, delim_whitespace=True, header=0, na_values=-998)
        positions = df['pos']
        Svec = []
        for i, row in df.iterrows():
            S = Stats(self.path2trained)
            for stat in S.stats:
                S.set_stat(stat, float(row[stat]))
            Svec.append(S)


        return [positions, Svec]


####### AODE emissions model #####################

    def _initialize_AODE(self,path2trained):
        if path2trained != '' and path2trained[-1] != '/':
            path2trained += '/'
        self.path2trained = path2trained

        file = open(self.path2trained+'component_stats.txt','r')
        f = file.read()
        file.close()
        f = f.strip().splitlines()
        self.statlist = [x.strip() for x in f]

        self.num2stat = {i:self.statlist[i] for i in range(len(self.statlist))}
        self.stat2num = {y:x for x,y in self.num2stat.items()}

        self.path2AODE = self.path2trained+'AODE_params/'   

        file = open(self.path2trained+'classes.txt','r')
        f = file.read()
        file.close()
        self.scenarios = [x.strip() for x in f.strip().splitlines()]
        print(self.scenarios)        

        self.JOINTS = [[[] for stat2 in self.stat2num.keys()] for stat1 in self.stat2num.keys()]
        self.MARGINALS = [[] for stat in self.stat2num.keys()]

        for stat1 in self.statlist:
            for stat2 in [x for x in self.statlist if x!= stat1]:
                statnum1 = self.stat2num[stat1]
                statnum2 = self.stat2num[stat2]
                if statnum1 < statnum2:
                    for scenario in self.scenarios:
                        self.JOINTS[statnum1][statnum2].append(self.gmm_fit(stat1,stat2,scenario))

        for stat1 in self.statlist:
            statnum1 = self.stat2num[stat1]
            for scenario in self.scenarios:
                self.MARGINALS[statnum1].append(self.gmm_fit_1D(stat1,scenario))

    def ode_likelihood(self,keystat,stats,scenario,imputation_mode='ode_compensation'): #adapted from ode function in SWIFr.py to only return likelihood for a single classification scenario
        score = stats.stat2score[keystat]
        print(score)
        if np.isnan(score):
            return 'n/a'
        else:
            Likelihood_list = []
            MARGS = self.MARGINALS[self.stat2num[keystat]]
            scenarionum = self.scenarios.index(scenario)
            M = MARGS[scenarionum]
            value = self.GMM_pdf(M,score)
            Likelihood = value
            Likelihood_list.append(value)
            stats_undefined = 0 #keep track of how many comparison stats are undefined
            for stat in self.statlist:
                if stat != keystat:
                    score2 = stats.stat2score[stat]
                    if np.isnan(score2):
                        stats_undefined += 1
                    else:
                        if self.stat2num[keystat] < self.stat2num[stat]:
                            H = self.conditional_GMM(score,2,self.JOINTS[self.stat2num[keystat]][self.stat2num[stat]][scenarionum])
                            value = self.GMM_pdf(H,score2)
                            Likelihood = Likelihood*value
                            Likelihood_list.append(value)
                        else:
                            H = self.conditional_GMM(score,1,self.JOINTS[self.stat2num[stat]][self.stat2num[keystat]][scenarionum])
                            value = self.GMM_pdf(H,score2)
                            Likelihood = Likelihood*value
                            Likelihood_list.append(value)
            if imputation_mode == 'ode_compensation':
                #print Likelihood_list
                n = stats_undefined
                xbar = float(sum(Likelihood_list))/len(Likelihood_list) #average conditional density for defined comparison stats
                for i in range(n):
                    Likelihood = Likelihood*xbar
            return Likelihood

    def GMM_pdf(self,G,x):
        w = G.weights_
        mu = G.means_
        C = G.covars_
        pdf = 0
        for i in range(len(w)):
            pdf += w[i]*self.normpdf(x,mu[i][0],math.sqrt(C[i][0]))
        return pdf


    def normpdf(self,x,mu,sigma):
        u = float(x-mu)/sigma
        y = (1/(math.sqrt(2*math.pi)*abs(sigma)))*math.exp(-u*u/2)
        return y
    
    def conditional_GMM(self,condval,keystat,G):
        #keystat = 1 if want stat1|stat2, keystat = 2 if want stat2|stat1
        H = Mix1D()

        for i in range(len(G.weights_)):
            sigma1 = math.sqrt(G.covars_[i][0][0])
            sigma2 = math.sqrt(G.covars_[i][1][1])
            ro = float(G.covars_[i][0][1])/(sigma1*sigma2)
            mu1 = G.means_[i][0]
            mu2 = G.means_[i][1]
            if keystat == 1:
                H.weights_.append(G.weights_[i])
                H.means_.append([mu1 + float(sigma1*ro*(condval-mu2))/sigma2])
                H.covars_.append([(1-ro**2)*sigma1**2])

            elif keystat == 2:
                H.weights_.append(G.weights_[i])
                H.means_.append([mu2 + float(sigma2*ro*(condval-mu1))/sigma1])
                H.covars_.append([(1-ro**2)*sigma2**2])
        return H

    def gmm_fit(self,stat1,stat2,scenario):
        G = pickle.load(open(self.path2AODE+stat1+'_'+stat2+'_'+scenario+'_GMMparams.p','rb'))
        return G

    def gmm_fit_1D(self,stat,scenario):
        G = pickle.load(open(self.path2AODE+stat+'_'+scenario+'_1D_GMMparams.p','rb'))
        return G

