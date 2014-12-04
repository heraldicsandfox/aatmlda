/*
 * Copyright (C) 2007 by
 * 
 * 	Xuan-Hieu Phan
 *	hieuxuan@ecei.tohoku.ac.jp or pxhieu@gmail.com
 * 	Graduate School of Information Sciences
 * 	Tohoku University
 *
 * GibbsLDA++ is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation; either version 2 of the License,
 * or (at your option) any later version.
 *
 * GibbsLDA++ is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GibbsLDA++; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
 */

/* 
 * References:
 * + The Java code of Gregor Heinrich (gregor@arbylon.net)
 *   http://www.arbylon.net/projects/LdaGibbsSampler.java
 * + "Parameter estimation for text analysis" by Gregor Heinrich
 *   http://www.arbylon.net/publications/text-est.pdf
 */

#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>
#include <time.h>
#include <float.h>
#include <set>
#include "constants.h"
#include "strtokenizer.h"
#include "utils.h"
#include "dataset.h"
#include "model.h"

using namespace std;

model::~model() {
    if (p) {
	delete p;
    }

    if (ptrndata) {
	delete ptrndata;
    }
    
    if (pnewdata) {
	delete pnewdata;
    }

    if (z) {
	for (int m = 0; m < M; m++) {
	    if (z[m]) {
		delete z[m];
	    }
	}
    }
    
    if (nw) {
	for (int w = 0; w < V; w++) {
	    if (nw[w]) {
		delete nw[w];
	    }
	}
    }

    if (word2constr) {
        delete word2constr;
    }

    if (constr2size) {
        delete constr2size;
    }

    if (nl) {
    for (int k = 0; k < K; k++) {
        if (nl[k]) {
        delete nl[k];
        }
    }
    }


    if (nd) {
	for (int m = 0; m < M; m++) {
	    if (nd[m]) {
		delete nd[m];
	    }
	}
    } 
    
    if (nwsum) {
	delete nwsum;
    }   
    
    if (ndsum) {
	delete ndsum;
    }
    
    if (theta) {
	for (int m = 0; m < M; m++) {
	    if (theta[m]) {
		delete theta[m];
	    }
	}
    }
    
    if (phi) {
	for (int k = 0; k < K; k++) {
	    if (phi[k]) {
		delete phi[k];
	    }
	}
    }

    // only for inference
    if (newz) {
	for (int m = 0; m < newM; m++) {
	    if (newz[m]) {
		delete newz[m];
	    }
	}
    }
    
    if (newnw) {
	for (int w = 0; w < newV; w++) {
	    if (newnw[w]) {
		delete newnw[w];
	    }
	}
    }

    if (newnd) {
	for (int m = 0; m < newM; m++) {
	    if (newnd[m]) {
		delete newnd[m];
	    }
	}
    } 
    
    if (newnwsum) {
	delete newnwsum;
    }   
    
    if (newndsum) {
	delete newndsum;
    }
    
    if (newtheta) {
	for (int m = 0; m < newM; m++) {
	    if (newtheta[m]) {
		delete newtheta[m];
	    }
	}
    }
    
    if (newphi) {
	for (int k = 0; k < K; k++) {
	    if (newphi[k]) {
		delete newphi[k];
	    }
	}
    }
}

void model::set_default_values() {
    wordmapfile = "wordmap.txt";
    trainlogfile = "trainlog.txt";
    tassign_suffix = ".tassign";
    theta_suffix = ".theta";
    phi_suffix = ".phi";
    others_suffix = ".others";
    twords_suffix = ".twords";
    
    dir = "./";
    dfile = "trndocs.dat";
    fname = "model";
    model_name = "model-final";    
    model_status = MODEL_STATUS_UNKNOWN;
    
    ptrndata = NULL;
    pnewdata = NULL;
    
    M = 0;
    V = 0;
    K = 100;
    alpha = 50.0 / K;
    beta = 0.1;
    niters = 2000;
    liter = 0;
    savestep = 200;    
    twords = 0;
    withrawstrs = 0;
    
    p = NULL;
    z = NULL;
    nw = NULL;
    nd = NULL;
    nl = NULL;
    nwsum = NULL;
    ndsum = NULL;
    theta = NULL;
    phi = NULL;
    word2constr = NULL;
    constr2size = NULL;
    
    newM = 0;
    newV = 0;
    newz = NULL;
    newnw = NULL;
    newnd = NULL;
    newnwsum = NULL;
    newndsum = NULL;
    newtheta = NULL;
    newphi = NULL;

    vacuous_dist = NULL;
    uniform_dist = NULL;
    doc_dists = NULL;
}

int model::parse_args(int argc, char ** argv) {
    return utils::parse_args(argc, argv, this);
}

int model::init(int argc, char ** argv) {
    // call parse_args
    if (parse_args(argc, argv)) {
	return 1;
    }
    
    if (model_status == MODEL_STATUS_EST) {
	// estimating the model from scratch
	if (init_est()) {
	    return 1;
	}
	
    } else if (model_status == MODEL_STATUS_ESTC) {
	// estimating the model from a previously estimated one
	if (init_estc()) {
	    return 1;
	}
	
    } else if (model_status == MODEL_STATUS_INF) {
	// do inference
	if (init_inf()) {
	    return 1;
	}
    }
    
    return 0;
}

int model::load_model(string model_name) {
    int i, j;
    
    string filename = dir + model_name + tassign_suffix;
    FILE * fin = fopen(filename.c_str(), "r");
    if (!fin) {
	printf("Cannot open file %s to load model!\n", filename.c_str());
	return 1;
    }
    
    char buff[BUFF_SIZE_LONG];
    string line;

    // allocate memory for z and ptrndata
    z = new int*[M];
    ptrndata = new dataset(M);
    ptrndata->V = V;

    for (i = 0; i < M; i++) {
	char * pointer = fgets(buff, BUFF_SIZE_LONG, fin);
	if (!pointer) {
	    printf("Invalid word-topic assignment file, check the number of docs!\n");
	    return 1;
	}
	
	line = buff;
	strtokenizer strtok(line, " \t\r\n");
	int length = strtok.count_tokens();
	
	vector<int> words;
	vector<int> topics;
	for (j = 0; j < length; j++) {
	    string token = strtok.token(j);
    
	    strtokenizer tok(token, ":");
	    if (tok.count_tokens() != 2) {
		printf("Invalid word-topic assignment line!\n");
		return 1;
	    }
	    
	    words.push_back(atoi(tok.token(0).c_str()));
	    topics.push_back(atoi(tok.token(1).c_str()));
	}
	
	// allocate and add new document to the corpus
	document * pdoc = new document(words);
	ptrndata->add_doc(pdoc, i);
	
	// assign values for z
	z[i] = new int[topics.size()];
	for (j = 0; j < topics.size(); j++) {
	    z[i][j] = topics[j];
	}
    }   
    
    fclose(fin);
    
    return 0;
}

int model::save_model(string model_name) {
    if (save_model_tassign(dir + model_name + tassign_suffix)) {
	return 1;
    }
    
    if (save_model_others(dir + model_name + others_suffix)) {
	return 1;
    }
    
    if (save_model_theta(dir + model_name + theta_suffix)) {
	return 1;
    }
    
    if (save_model_phi(dir + model_name + phi_suffix)) {
	return 1;
    }
    
    if (twords > 0) {
	if (save_model_twords(dir + model_name + twords_suffix)) {
	    return 1;
	}
    }
    
    return 0;
}

int model::save_model_tassign(string filename) {
    int i, j;
    
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }

    // wirte docs with topic assignments for words
    for (i = 0; i < ptrndata->M; i++) {    
	for (j = 0; j < ptrndata->docs[i]->length; j++) {
	    fprintf(fout, "%d:%d ", ptrndata->docs[i]->words[j], z[i][j]);
	}
	fprintf(fout, "\n");
    }

    fclose(fout);
    
    return 0;
}

int model::save_model_theta(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    for (int i = 0; i < M; i++) {
	for (int j = 0; j < K; j++) {
	    fprintf(fout, "%f ", theta[i][j]);
	}
	fprintf(fout, "\n");
    }
    
    fclose(fout);
    
    return 0;
}

int model::save_model_phi(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    for (int i = 0; i < K; i++) {
	for (int j = 0; j < V; j++) {
	    fprintf(fout, "%f ", phi[i][j]);
	}
	fprintf(fout, "\n");
    }
    
    fclose(fout);    
    
    return 0;
}

int model::save_model_others(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }

    fprintf(fout, "alpha=%f\n", alpha);
    fprintf(fout, "beta=%f\n", beta);
    fprintf(fout, "ntopics=%d\n", K);
    fprintf(fout, "ndocs=%d\n", M);
    fprintf(fout, "nwords=%d\n", V);
    fprintf(fout, "liter=%d\n", liter);
    
    fclose(fout);    
    
    return 0;
}

int model::save_model_twords(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    if (twords > V) {
	twords = V;
    }
    mapid2word::iterator it;
    
    for (int k = 0; k < K; k++) {
	vector<pair<int, double> > words_probs;
	pair<int, double> word_prob;
	for (int w = 0; w < V; w++) {
	    word_prob.first = w;
	    word_prob.second = phi[k][w];
	    words_probs.push_back(word_prob);
	}
    
        // quick sort to sort word-topic probability
	utils::quicksort(words_probs, 0, words_probs.size() - 1);
	
	fprintf(fout, "Topic %dth:\n", k);
	for (int i = 0; i < twords; i++) {
	    it = id2word.find(words_probs[i].first);
	    if (it != id2word.end()) {
		fprintf(fout, "\t%s   %f\n", (it->second).c_str(), words_probs[i].second);
	    }
	}
    }
    
    fclose(fout);    
    
    return 0;    
}

int model::save_inf_model(string model_name) {
    if (save_inf_model_tassign(dir + model_name + tassign_suffix)) {
	return 1;
    }
    
    if (save_inf_model_others(dir + model_name + others_suffix)) {
	return 1;
    }
    
    if (save_inf_model_newtheta(dir + model_name + theta_suffix)) {
	return 1;
    }
    
    if (save_inf_model_newphi(dir + model_name + phi_suffix)) {
	return 1;
    }

    if (twords > 0) {
	if (save_inf_model_twords(dir + model_name + twords_suffix)) {
	    return 1;
	}
    }
    
    return 0;
}

int model::save_inf_model_tassign(string filename) {
    int i, j;
    
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }

    // wirte docs with topic assignments for words
    for (i = 0; i < pnewdata->M; i++) {    
	for (j = 0; j < pnewdata->docs[i]->length; j++) {
	    fprintf(fout, "%d:%d ", pnewdata->docs[i]->words[j], newz[i][j]);
	}
	fprintf(fout, "\n");
    }

    fclose(fout);
    
    return 0;
}

int model::save_inf_model_newtheta(string filename) {
    int i, j;

    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    for (i = 0; i < newM; i++) {
	for (j = 0; j < K; j++) {
	    fprintf(fout, "%f ", newtheta[i][j]);
	}
	fprintf(fout, "\n");
    }
    
    fclose(fout);
    
    return 0;
}

int model::save_inf_model_newphi(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    for (int i = 0; i < K; i++) {
	for (int j = 0; j < newV; j++) {
	    fprintf(fout, "%f ", newphi[i][j]);
	}
	fprintf(fout, "\n");
    }
    
    fclose(fout);    
    
    return 0;
}

int model::save_inf_model_others(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }

    fprintf(fout, "alpha=%f\n", alpha);
    fprintf(fout, "beta=%f\n", beta);
    fprintf(fout, "eta=%f\n", eta);
    fprintf(fout, "ntopics=%d\n", K);
    fprintf(fout, "ndocs=%d\n", newM);
    fprintf(fout, "nwords=%d\n", newV);
    fprintf(fout, "liter=%d\n", inf_liter);
    
    fclose(fout);    
    
    return 0;
}

int model::save_inf_model_twords(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    if (twords > newV) {
	twords = newV;
    }
    mapid2word::iterator it;
    map<int, int>::iterator _it;
    
    for (int k = 0; k < K; k++) {
	vector<pair<int, double> > words_probs;
	pair<int, double> word_prob;
	for (int w = 0; w < newV; w++) {
	    word_prob.first = w;
	    word_prob.second = newphi[k][w];
	    words_probs.push_back(word_prob);
	}
    
        // quick sort to sort word-topic probability
	utils::quicksort(words_probs, 0, words_probs.size() - 1);
	
	fprintf(fout, "Topic %dth:\n", k);
	for (int i = 0; i < twords; i++) {
	    _it = pnewdata->_id2id.find(words_probs[i].first);
	    if (_it == pnewdata->_id2id.end()) {
		continue;
	    }
	    it = id2word.find(_it->second);
	    if (it != id2word.end()) {
		fprintf(fout, "\t%s   %f\n", (it->second).c_str(), words_probs[i].second);
	    }
	}
    }
    
    fclose(fout);    
    
    return 0;    
}

int model::save_model_words(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
    printf("Cannot open file %s to save!\n", filename.c_str());
    return 1;
    }

    mapid2word::iterator it = id2word.begin();
    while (it != id2word.end()) {
        fprintf(fout, "%d   %s\n", (it->first), (it->second).c_str());
        ++it;
    }

    fclose(fout);
    return 0;
}


int model::init_est() {
    int m, n, w, k;

    p = new double[K];

    // + read training data
    ptrndata = new dataset;
    if (ptrndata->read_trndata(dir + dfile, dir + wordmapfile, dir + confile)) {
        printf("Fail to read training data!\n");
        return 1;
    }
		
    // + allocate memory and assign values for variables
    M = ptrndata->M;
    V = ptrndata->V;
    C = ptrndata->C;
    // K: from command line or default value
    // alpha, beta, eta: from command line or default values
    // niters, savestep: from command line or default values

    nw = new int*[V];
    for (w = 0; w < V; w++) {
        nw[w] = new int[K];
        for (k = 0; k < K; k++) {
    	    nw[w][k] = 0;
        }
    }
	
    nd = new int*[M];
    for (m = 0; m < M; m++) {
        nd[m] = new int[K];
        for (k = 0; k < K; k++) {
    	    nd[m][k] = 0;
        }
    }

    nl = new int*[K];
    for (k = 0; k < K; k++) {
        nl[k] = new int[C];
        for (int l = 0; l < C; ++l) {
            nl[k][l] = 0;
        }
    }
	
    nwsum = new int[K];
    for (k = 0; k < K; k++) {
	nwsum[k] = 0;
    }
    
    ndsum = new int[M];
    for (m = 0; m < M; m++) {
	ndsum[m] = 0;
    }

    constr2size = new int[C];
    for (int l = 0; l < C; ++l) {
        constr2size[l] = 0;
    }

    word2constr = new int[V];
    for (w = 0; w < V; ++w) {
        map<int, int>::iterator it = ptrndata->word2constraint.find(w);
        if (it != ptrndata->word2constraint.end()) {
            word2constr[w] = it->second;
            constr2size[it->second]++;
        } else {
            word2constr[w] = -1;
        }
    }

    doc_dists = new double*[M];
    for (m = 0; m < M; m++) {
        doc_dists[m] = new double[V];
    }
    
    vacuous_dist = new double[V];
    uniform_dist = new double[V];
    
    for (w = 0; w < V; ++w) {
        vacuous_dist[w] = 0.0;
        uniform_dist[w] = 1.0 / V;
        for (m = 0; m < M; ++m) {
            doc_dists[m][w] = 0.0;
        }
    }

    int numtokens = 0;
    srandom(time(0)); // initialize for random number generation
    z = new int*[M];
    for (m = 0; m < ptrndata->M; m++) {
	int N = ptrndata->docs[m]->length;
	z[m] = new int[N];
    map<int, int>::iterator it;
    
        // initialize for z
        for (n = 0; n < N; n++) {
    	    int topic = (int)(((double)random() / RAND_MAX) * K);
    	    z[m][n] = topic;
    	    
            int tmpwd = ptrndata->docs[m]->words[n];

    	    // number of instances of word i assigned to topic j
    	    nw[tmpwd][topic] += 1;
    	    // number of words in document i assigned to topic j
    	    nd[m][topic] += 1;
    	    // total number of words assigned to topic j
    	    nwsum[topic] += 1;

            if (word2constr[tmpwd] > -1) {
                nl[topic][word2constr[tmpwd]]++;
            }


            vacuous_dist[tmpwd] += 1.0;
            doc_dists[m][tmpwd] += 1.0;
        } 
        // total number of words in document i
        ndsum[m] = N; 
        numtokens += N;     
    }
    
    // This gives us final distributions for vacuous, uniform, and doc-y
    for (w = 0; w < V; ++w) {
        vacuous_dist[w] = (vacuous_dist[w] + beta) / (numtokens + V * beta);
        for (m = 0; m < M; ++m) {
            doc_dists[m][w] = (doc_dists[m][w] + beta) / (numtokens + V * beta);
        }
    }

	
    phi = new double*[K];
    for (k = 0; k < K; k++) {
        phi[k] = new double[V];
    }    

    theta = new double*[M];
    for (m = 0; m < M; m++) {
        theta[m] = new double[K];
    }
    
    return 0;
}

int model::init_estc() {
    // estimating the model from a previously estimated one
    int m, n, w, k;

    p = new double[K];

    // load moel, i.e., read z and ptrndata
    if (load_model(model_name)) {
	printf("Fail to load word-topic assignmetn file of the model!\n");
	return 1;
    }

    nw = new int*[V];
    for (w = 0; w < V; w++) {
        nw[w] = new int[K];
        for (k = 0; k < K; k++) {
    	    nw[w][k] = 0;
        }
    }
	
    nd = new int*[M];
    for (m = 0; m < M; m++) {
        nd[m] = new int[K];
        for (k = 0; k < K; k++) {
    	    nd[m][k] = 0;
        }
    }
	
    nwsum = new int[K];
    for (k = 0; k < K; k++) {
	nwsum[k] = 0;
    }
    
    ndsum = new int[M];
    for (m = 0; m < M; m++) {
	ndsum[m] = 0;
    }

    for (m = 0; m < ptrndata->M; m++) {
	int N = ptrndata->docs[m]->length;

	// assign values for nw, nd, nwsum, and ndsum
        for (n = 0; n < N; n++) {
    	    int w = ptrndata->docs[m]->words[n];
    	    int topic = z[m][n];
    	    
    	    // number of instances of word i assigned to topic j
    	    nw[w][topic] += 1;
    	    // number of words in document i assigned to topic j
    	    nd[m][topic] += 1;
    	    // total number of words assigned to topic j
    	    nwsum[topic] += 1;
        } 
        // total number of words in document i
        ndsum[m] = N;
    }
	
    theta = new double*[M];
    for (m = 0; m < M; m++) {
        theta[m] = new double[K];
    }
	
    phi = new double*[K];
    for (k = 0; k < K; k++) {
        phi[k] = new double[V];
    }    

    return 0;
}

void model::estimate() {
    if (twords > 0) {
	// print out top words per topic
	dataset::read_wordmap(dir + wordmapfile, &id2word);
    }

    if (save_model_words(dir + fname + ".words")) {
        printf("Could not create word file!");
        return;
    }

    printf("Sampling %d iterations!\n", niters);

    int last_iter = liter;
    for (liter = last_iter + 1; liter <= niters + last_iter; liter++) {
        printf("Iteration %d ...\n", liter);
    	
    	// for all z_i (Each Gibbs sample)
        for (int m = 0; m < M; m++) {
            for (int n = 0; n < ptrndata->docs[m]->length; n++) {
    	        // (z_i = z[m][n])
    	        // sample from p(z_i|z_-i, w)
    	        int topic = sampling(m, n);
    	        z[m][n] = topic;
    	    }
    	}
    	
    	if (savestep > 0) {
    	    if (liter % savestep == 0 && liter != niters + last_iter) {
    		    // saving the model
    		    printf("Saving the model at iteration %d ...\n", liter);
    		    compute_theta();
    		    compute_phi();
    		    save_model(utils::generate_model_name(liter, fname));
                printf("Starting AATM...\n");
                aatm();
                printf("Saving post-AATM...\n");
                compute_theta();
                compute_phi();
                save_model(utils::generate_model_name_aatm(liter, fname));
    	    }
    	}
    }
    printf("Gibbs sampling completed!\n");
    printf("Saving the final model!\n");
    compute_theta();
    compute_phi();
    liter--;
    save_model(utils::generate_model_name(-1, fname));
}

void model::aatm() {
    // This will keep track of deleted topics
    set<int> freetopics;

    // Figure out which topics should be dumped using KL divergence
    double * vdivergences = new double[K];
    double top_vdivergence = 0.0;
    double * udivergences = new double[K];
    double top_udivergence = 0.0;

    for (int k = 0; k < K; ++k) {
        vdivergences[k] = 0.0;
        for (int w = 0; w < V; ++w) {
            vdivergences[k] += (phi[k][w] - vacuous_dist[w]) * log(phi[k][w] / vacuous_dist[w]);
            udivergences[k] += (phi[k][w] - uniform_dist[w]) * log(phi[k][w] / uniform_dist[w]);
        }
        if (vdivergences[k] > top_vdivergence) {
            top_vdivergence = vdivergences[k];
        }
        if (udivergences[k] > top_udivergence) {
            top_udivergence = udivergences[k];
        }
    }

    double u_threshold = 0.1 * top_udivergence;
    double v_threshold = 0.1 * top_vdivergence;

    for (int k = 0; k < K; ++k) {
        if (udivergences[k] < u_threshold || vdivergences[k] < v_threshold) {
            delete_topic(k);
            freetopics.insert(k);
            printf("Deleting topic %d ...\n", k);
        }
    }

    // Figure out which topics are similar enough to join
    double ** kdivergences = new double*[K];
    double top_kdivergence = 0.0;
    for (int k = 0; k < K; ++k) {
        kdivergences[k] = new double[K];
        if (freetopics.find(k) != freetopics.end())
            continue;
        for (int l = k + 1; l < K; ++l) {
            if (freetopics.find(l) != freetopics.end())
                continue;
            kdivergences[k][l] = 0.0;
            for (int w = 0; w < V; ++w) {
                kdivergences[k][l] += (phi[k][w] - phi[l][w]) * log(phi[k][w] / phi[l][w]);
            }
            if (kdivergences[k][l] > top_kdivergence) {
                top_kdivergence = kdivergences[k][l];
            }
        }
    }

    double k_threshold = 0.05 * top_kdivergence;

    for (int k = 0; k < K; ++k) {
        if (freetopics.find(k) != freetopics.end())
            continue;
        for (int l = k + 1; l < K; ++l) {
            if (freetopics.find(k) != freetopics.end())
                continue;
            if (kdivergences[k][l] < k_threshold && freetopics.size() < ((K + 1)/2)) {
                merge_topics(k, l);
                freetopics.insert(l);
                printf("Deleting topic %d ...\n", k);
            }
        }
    }

    int topics_to_create = freetopics.size();

    // Figure out which topics are doc-y enough to be split
    double ** mdivergences = new double*[K];
    double top_mdivergence = 0.0;
    for (int k = 0; k < K; ++k) {
        mdivergences[k] = new double[M];
        for (int m = 0; m < M; ++m) {
            if (freetopics.find(k) != freetopics.end()) {
                mdivergences[k][m] = DBL_MAX;
            } else {
                mdivergences[k][m] = 0.0;
                for (int w = 0; w < V; ++w) {
                    mdivergences[k][m] += (phi[k][w] - doc_dists[m][w]) * log(phi[k][w] / doc_dists[m][w]);
                }
            }
        }
    }
    vector<int> topics_to_split = utils::split_topics(kdivergences, K, M, topics_to_create);
    set<int>::iterator ft = freetopics.begin();
    vector<int>::iterator ts = topics_to_split.begin();
    while (ft != freetopics.end() && ts != topics_to_split.end()) {
        printf("Splitting topic %d with %d ...\n", *ts, *ft);
        split_topic(*ft, *ts);
        ft++;
        ts++;
    }

    delete [] vdivergences;
    delete [] udivergences;
    for (int k = 0; k < K; ++k) {
        delete [] kdivergences[k];
        delete [] mdivergences[k];
    }
    delete [] kdivergences;
    delete [] mdivergences;
}

void model::delete_topic(int top) {
    for (int m = 0; m < M; m++) {
        for (int n = 0; n < ptrndata->docs[m]->length; n++) {
            if (z[m][n] == top) {
                int w = ptrndata->docs[m]->words[n];
                z[m][n] = -1;
                nw[w][top] -= 1;
                nd[m][top] -= 1;
                nwsum[top] -= 1;
                ndsum[m] -= 1;
            }
        }
    }
}

void model::split_topic(int freet, int splitt) {
    for (int m = 0; m < M; m++) {
        for (int n = 0; n < ptrndata->docs[m]->length; n++) {
            if (z[m][n] == splitt && ((double)random() / RAND_MAX) > 0.5) {
                int w = ptrndata->docs[m]->words[n];
                z[m][n] = freet;
                nw[w][splitt] -= 1;
                nd[m][splitt] -= 1;
                nwsum[splitt] -= 1;
                nw[w][freet] += 1;
                nd[m][freet] += 1;
                nwsum[freet] += 1;
            }
        }
    }
}

void model::merge_topics(int newt, int oldt) {
    for (int m = 0; m < M; m++) {
        for (int n = 0; n < ptrndata->docs[m]->length; n++) {
            if (z[m][n] == oldt) {
                int w = ptrndata->docs[m]->words[n];
                z[m][n] = newt;
                nw[w][oldt] -= 1;
                nd[m][oldt] -= 1;
                nwsum[oldt] -= 1;
                nw[w][newt] += 1;
                nd[m][newt] += 1;
                nwsum[newt] += 1;
            }
        }
    }
}

// The Gibbs sampler itself
int model::sampling(int m, int n) {
    // remove z_i from the count variables
    int topic = z[m][n];
    int w = ptrndata->docs[m]->words[n];
    int l = word2constr[w];
    if (topic != -1) {
        nw[w][topic] -= 1;
        nd[m][topic] -= 1;
        nwsum[topic] -= 1;
        ndsum[m] -= 1;
        if (l != -1) {
            nl[topic][l] -= 1;
        }
    }
    
    double Vbeta = V * beta;
    double Kalpha = K * alpha;
    if (l == -1) {  
        // do multinomial sampling via cumulative method
        for (int k = 0; k < K; k++) {
	        p[k] = (nw[w][k] + beta) / (nwsum[k] + Vbeta) *
		        (nd[m][k] + alpha) / (ndsum[m] + Kalpha);
        }
        // cumulate multinomial parameters
        for (int k = 1; k < K; k++) {
	        p[k] += p[k - 1];
        }
        // scaled sample because of unnormalized p[]
        double u = ((double)random() / RAND_MAX) * p[K - 1];
        
        for (topic = 0; topic < K; topic++) {
    	if (p[topic] > u) {
    	    break;
    	}
        }
    } else {
        // do multinomial sampling via cumulative method
        int csize = constr2size[l];
        double Cbeta = csize * beta;
        double Ceta = csize * eta;
        for (int k = 0; k < K; k++) {
            p[k] = (nd[m][k] + alpha) / (ndsum[m] + Kalpha) *
                   (nl[k][l] + Cbeta) / (nwsum[k] + Vbeta) *
                   (nw[w][k] + eta) / (nl[k][l] + Ceta); 
        }
        // cumulate multinomial parameters
        for (int k = 1; k < K; k++) {
            p[k] += p[k - 1];
        }
        // scaled sample because of unnormalized p[]
        double u = ((double)random() / RAND_MAX) * p[K - 1];
        
        for (topic = 0; topic < K; topic++) {
        if (p[topic] > u) {
            break;
        }
        }
        // count the special constraint value
        nl[topic][l] += 1;
    }
    
    // add newly estimated z_i to count variables
    nw[w][topic] += 1;
    nd[m][topic] += 1;
    nwsum[topic] += 1;
    ndsum[m] += 1;    
    
    return topic;
}

void model::compute_theta() {
    for (int m = 0; m < M; m++) {
	for (int k = 0; k < K; k++) {
	    theta[m][k] = (nd[m][k] + alpha) / (ndsum[m] + K * alpha);
	}
    }
}

void model::compute_phi() {
    for (int k = 0; k < K; k++) {
	for (int w = 0; w < V; w++) {
        if (word2constr[w] == -1) {
	        phi[k][w] = (nw[w][k] + beta) / (nwsum[k] + V * beta);
        } else {
            int l = word2constr[w];
            int csize = constr2size[l];
            phi[k][w] = (nl[k][l] + csize * beta) / (nwsum[k] + V * beta) *
                        (nw[w][k] + eta) / (nl[k][l] + csize * eta); 
        }
	}
    }
}

int model::init_inf() {
    // estimating the model from a previously estimated one
    int m, n, w, k;

    p = new double[K];

    // load moel, i.e., read z and ptrndata
    if (load_model(model_name)) {
	printf("Fail to load word-topic assignmetn file of the model!\n");
	return 1;
    }

    nw = new int*[V];
    for (w = 0; w < V; w++) {
        nw[w] = new int[K];
        for (k = 0; k < K; k++) {
    	    nw[w][k] = 0;
        }
    }
	
    nd = new int*[M];
    for (m = 0; m < M; m++) {
        nd[m] = new int[K];
        for (k = 0; k < K; k++) {
    	    nd[m][k] = 0;
        }
    }
	
    nwsum = new int[K];
    for (k = 0; k < K; k++) {
	nwsum[k] = 0;
    }
    
    ndsum = new int[M];
    for (m = 0; m < M; m++) {
	ndsum[m] = 0;
    }

    for (m = 0; m < ptrndata->M; m++) {
	int N = ptrndata->docs[m]->length;

	// assign values for nw, nd, nwsum, and ndsum	
        for (n = 0; n < N; n++) {
    	    int w = ptrndata->docs[m]->words[n];
    	    int topic = z[m][n];
    	    
    	    // number of instances of word i assigned to topic j
    	    nw[w][topic] += 1;
    	    // number of words in document i assigned to topic j
    	    nd[m][topic] += 1;
    	    // total number of words assigned to topic j
    	    nwsum[topic] += 1;
        } 
        // total number of words in document i
        ndsum[m] = N;      
    }
    
    // read new data for inference
    pnewdata = new dataset;
    if (withrawstrs) {
	if (pnewdata->read_newdata_withrawstrs(dir + dfile, dir + wordmapfile)) {
    	    printf("Fail to read new data!\n");
    	    return 1;
	}    
    } else {
	if (pnewdata->read_newdata(dir + dfile, dir + wordmapfile)) {
    	    printf("Fail to read new data!\n");
    	    return 1;
	}    
    }
    
    newM = pnewdata->M;
    newV = pnewdata->V;
    
    newnw = new int*[newV];
    for (w = 0; w < newV; w++) {
        newnw[w] = new int[K];
        for (k = 0; k < K; k++) {
    	    newnw[w][k] = 0;
        }
    }
	
    newnd = new int*[newM];
    for (m = 0; m < newM; m++) {
        newnd[m] = new int[K];
        for (k = 0; k < K; k++) {
    	    newnd[m][k] = 0;
        }
    }
	
    newnwsum = new int[K];
    for (k = 0; k < K; k++) {
	newnwsum[k] = 0;
    }
    
    newndsum = new int[newM];
    for (m = 0; m < newM; m++) {
	newndsum[m] = 0;
    }

    srandom(time(0)); // initialize for random number generation
    newz = new int*[newM];
    for (m = 0; m < pnewdata->M; m++) {
	int N = pnewdata->docs[m]->length;
	newz[m] = new int[N];

	// assign values for nw, nd, nwsum, and ndsum	
        for (n = 0; n < N; n++) {
    	    int w = pnewdata->docs[m]->words[n];
    	    int _w = pnewdata->_docs[m]->words[n];
    	    int topic = (int)(((double)random() / RAND_MAX) * K);
    	    newz[m][n] = topic;
    	    
    	    // number of instances of word i assigned to topic j
    	    newnw[_w][topic] += 1;
    	    // number of words in document i assigned to topic j
    	    newnd[m][topic] += 1;
    	    // total number of words assigned to topic j
    	    newnwsum[topic] += 1;
        } 
        // total number of words in document i
        newndsum[m] = N;      
    }    
    
    newtheta = new double*[newM];
    for (m = 0; m < newM; m++) {
        newtheta[m] = new double[K];
    }
	
    newphi = new double*[K];
    for (k = 0; k < K; k++) {
        newphi[k] = new double[newV];
    }    
    
    return 0;        
}

void model::inference() {
    if (twords > 0) {
	// print out top words per topic
	dataset::read_wordmap(dir + wordmapfile, &id2word);
    }

    printf("Sampling %d iterations for inference!\n", niters);
    
    for (inf_liter = 1; inf_liter <= niters; inf_liter++) {
	printf("Iteration %d ...\n", inf_liter);
	
	// for all newz_i
	for (int m = 0; m < newM; m++) {
	    for (int n = 0; n < pnewdata->docs[m]->length; n++) {
		// (newz_i = newz[m][n])
		// sample from p(z_i|z_-i, w)
		int topic = inf_sampling(m, n);
		newz[m][n] = topic;
	    }
	}
    }
    
    printf("Gibbs sampling for inference completed!\n");
    printf("Saving the inference outputs!\n");
    compute_newtheta();
    compute_newphi();
    inf_liter--;
    save_inf_model(dfile);
}

int model::inf_sampling(int m, int n) {
    // remove z_i from the count variables
    int topic = newz[m][n];
    int w = pnewdata->docs[m]->words[n];
    int _w = pnewdata->_docs[m]->words[n];
    newnw[_w][topic] -= 1;
    newnd[m][topic] -= 1;
    newnwsum[topic] -= 1;
    newndsum[m] -= 1;
    
    double Vbeta = V * beta;
    double Kalpha = K * alpha;
    // do multinomial sampling via cumulative method
    for (int k = 0; k < K; k++) {
	p[k] = (nw[w][k] + newnw[_w][k] + beta) / (nwsum[k] + newnwsum[k] + Vbeta) *
		    (newnd[m][k] + alpha) / (newndsum[m] + Kalpha);
    }
    // cumulate multinomial parameters
    for (int k = 1; k < K; k++) {
	p[k] += p[k - 1];
    }
    // scaled sample because of unnormalized p[]
    double u = ((double)random() / RAND_MAX) * p[K - 1];
    
    for (topic = 0; topic < K; topic++) {
	if (p[topic] > u) {
	    break;
	}
    }
    
    // add newly estimated z_i to count variables
    newnw[_w][topic] += 1;
    newnd[m][topic] += 1;
    newnwsum[topic] += 1;
    newndsum[m] += 1;    
    
    return topic;
}

void model::compute_newtheta() {
    for (int m = 0; m < newM; m++) {
	for (int k = 0; k < K; k++) {
	    newtheta[m][k] = (newnd[m][k] + alpha) / (newndsum[m] + K * alpha);
	}
    }
}

void model::compute_newphi() {
    map<int, int>::iterator it;
    for (int k = 0; k < K; k++) {
	for (int w = 0; w < newV; w++) {
	    it = pnewdata->_id2id.find(w);
	    if (it != pnewdata->_id2id.end()) {
		newphi[k][w] = (nw[it->second][k] + newnw[w][k] + beta) / (nwsum[k] + newnwsum[k] + V * beta);
	    }
	}
    }
}

