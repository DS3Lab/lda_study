# lda_study

study various inference methods for latent dirichlet allocation, including collapsed gibbs sampling, variational inference, considering both single and distributed implementation.

### Gibbs Sampler
[FastLDA: Fast Collapsed Gibbs Sampling For Latent Dirichlet
Allocation, KDD'2008](http://www.ics.uci.edu/~asuncion/pubs/KDD_08.pdf)

[SparseLDA: Efficient Methods for Topic Model Inference on Streaming
Document Collections, KDD'2009](https://core.ac.uk/download/pdf/21747811.pdf)

[AliasLDA: Reducing the Sampling Complexity of Topic Models, KDD'2014](https://pdfs.semanticscholar.org/137a/ec8c56102cea1ac7c083989036bb51331fdc.pdf)

[F+LDA: A Scalable Asynchronous Distributed Algorithm for Topic Modeling, 2014](http://arxiv.org/pdf/1412.4986v1.pdf)

[LightLDA: Big Topic Models on Modest Compute Clusters, 2014](https://arxiv.org/pdf/1412.1576v1.pdf)

[WarpLDA: a Cache Efficient O(1) Algorithm for Latent Dirichlet Allocation, VLDB'2016](http://www.vldb.org/pvldb/vol9/p744-chen.pdf)

### Variational inference

[CVB: A Collapsed Variational Bayesian Inference Algorithm for Latent Dirichlet Allocation, 2007](http://papers.nips.cc/paper/3113-a-collapsed-variational-bayesian-inference-algorithm-for-latent-dirichlet-allocation.pdf)

[Stochastic Variational Inference, 2013](http://www.jmlr.org/papers/volume14/hoffman13a/hoffman13a.pdf)
	
	Stochastic method can avoid scanning the whole dataset at each iteration, which is time-consuming in batch mode.
	It iterates between subsampling of data and adjusting the hidden structure based only on the subsample.
	Variational inference is amenable to stochastic optimization because the variational objective decomposes into a sum of terms, one for each data point in the analysis.

[Stochastic Collapsed Variational Bayesian Inference for
Latent Dirichlet Allocation, KDD'2013](http://www.ics.uci.edu/~welling/publications/papers/fp1199-foulds.pdf)

[CVB0: On Smoothing and Inference for Topic Models, 2009](http://arxiv.org/pdf/1205.2662.pdf)

[Sparse stochastic inference for latent Dirichlet allocation, ICML'2012](http://www.cs.columbia.edu/~blei/papers/MimnoHoffmanBlei2012.pdf)

### Stochatic gradient Sampler
[Distributing the Stochastic Gradient Sampler for
Large-Scale LDA, KDD'2016](http://www.kdd.org/kdd2016/papers/files/rpp0277-yangA.pdf)

### Belief propagation
[Learning Topic Models by Belief Propagation, 2012](https://arxiv.org/pdf/1109.3437.pdf)

[Residual Belief Propagation for Topic Modeling, 2012](http://arxiv.org/pdf/1204.6610.pdf)

### Others
[Memory-Efficient Topic Modeling](http://arxiv.org/pdf/1206.1147.pdf)

### Distributed Implementaion

#### variational method
[Mr. LDA: A Flexible Large Scale Topic Modeling Package
using Variational Inference in MapReduce, WWW'2012](http://kzhai.github.io/paper/2012_www.pdf)

#### gibbs sampling
Two types of implementation

- Share Word-Topic matirx using PS
	
	[YahooLDA: An Architecture for Parallel Topic Models](http://vldb.org/pvldb/vldb2010/papers/R63.pdf)
	
- Shuffle Doc-Word matrix with topic assignment 
	
	[WarpLDA: a Cache Efficient O(1) Algorithm for Latent Dirichlet Allocation, VLDB'2016](http://www.vldb.org/pvldb/vol9/p744-chen.pdf)

### Study papers
[On Smoothing and Inference for Topic Models, 2009](http://arxiv.org/pdf/1205.2662.pdf)

[Rethinking Collapsed Variational Bayes Inference for LDA, ICML'2012](http://icml.cc/2012/papers/530.pdf)

[Variational Inference: A Review for Statisticians](https://arxiv.org/pdf/1601.00670v3.pdf)


### Some discuss
[Variational vs MCMC, quora](https://www.quora.com/When-should-I-prefer-variational-inference-over-MCMC-for-Bayesian-analysis)

### Reading List
[LDA Reading List](http://bigml.cs.tsinghua.edu.cn/~jianfei/lda-reading.html)