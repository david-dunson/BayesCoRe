require("msm")
require("gtools")

options(max.print =1000)
GibbsSampler<- function(steps){

	mu1 = c(d:1)
	mu2 = c(d:1)

	rmul<- function(p){
		c(1:length(p))[rmultinom(1,1,p)==1]
	}

	sigma1 = 0.1
	sigma2 = 0.1

	p = rep(1/d,d)

	trace_p = numeric()
	# trace_sigma =numeric()
	trace_mu1 = numeric()
	trace_mu2 = numeric()
	trace_Z = numeric()

	for(j in 1:steps){
		#update latent Z
		diff1 = outer(y1,mu1,"-")
		diff2 = outer(y2,mu2,"-")

		normal_density = - diff1^2 / 2/sigma1^2 - diff2^2 / 2/sigma2^2

		p_z = matrix(log(p),N,d,byrow=T) + normal_density

		p_z = t(apply(p_z, 1, function(x){
			x-max(x)
			}))

		Z = apply(exp(p_z), 1,rmul)

		Z_count = numeric(d)
		#update mu

		for(i in 1:d){
			
			# if(i==1){
			# 	# ub = 100
			# 	# }else{
			# 	# 	ub = mu[i-1]
			# 	# }

			# 	# if(i==d){
			# 	# 	lb = 0
			# 	# 	}else{
			# 	# 		lb = mu[i+1]
			# 	# 	}

					n_sum = sum(Z==i)
					Z_count[i] = n_sum

					y1_sum = sum(y1[Z==i])
					v = 1/(n_sum/sigma1^2 + 1/10)
					m = v* y1_sum/sigma1^2
					mu1[i] = rnorm(1,m,sqrt(v))

					y2_sum = sum(y2[Z==i])

					v = 1/(n_sum/sigma2^2 + 1/10)
					m = v* y2_sum/sigma2^2
					mu2[i] = rnorm(1,m,sqrt(v))


				}


		#update sigma

		a = N/2 + 2

		diff = y1- mu1[Z]
		b = sum(diff^2)/2 + 1

		sigma1 = sqrt(1/ rgamma(1,a,b))


		diff = y2- mu2[Z]
		b = sum(diff^2)/2 + 1

		sigma2 = sqrt(1/ rgamma(1,a,b))

		#update p
		# alpha = rep(1, d)
		alpha1 = rep(alpha,d)

		for(i in 1:d){
			alpha1[i] = alpha + sum(Z==i)
		}

		p = c(rdirichlet(1, alpha1))

		trace_p = rbind(trace_p,p)
		trace_mu1 = rbind(trace_mu1,mu1)
		trace_mu2 = rbind(trace_mu2,mu2)
		trace_Z = rbind(trace_Z, Z_count)
		# trace_sigma = c(trace_sigma, sigma)	

		print(j)
}

	list("p" = trace_p,"mu1"= trace_mu1,"mu2"= trace_mu2, "Z" = trace_Z)
}


d = 3

alpha =2

Gibbsfit= GibbsSampler(10000)

sampling_idx_gibbs<- c(5001:10000)

require("ggplot2")

order= c(2,3,1)

# save(Gibbsfit, file="Gibbsfit.RDa")


pdf("../draft/fmm_mu2_gibbs.pdf",6,3)
df = data.frame( step =c(1:5000), mu2=c(Gibbsfit$mu2[sampling_idx_gibbs,order]), component =as.factor(rep(c(1:3), each=5000)))
ggplot(data=df) + geom_path( aes(x=step, y=mu2, col=component),alpha =0.8,size=.8) + theme_bw()+ylim(0,6)
dev.off()



pdf("../draft/fmm_mu1_gibbs.pdf",6,3)
df = data.frame( step =c(1:5000), mu1=c(Gibbsfit$mu1[sampling_idx_gibbs,order]), component =as.factor(rep(c(1:3), each=5000)))
ggplot(data=df) + geom_path( aes(x=step, y=mu1, col=component),alpha =0.8,size=.8) + theme_bw()+ylim(0,6)
dev.off()



pdf("../draft/fmm_w_gibbs.pdf",6,3)
df = data.frame( step =c(1:5000), w=c(Gibbsfit$p[sampling_idx_gibbs,order]), component =as.factor(rep(c(1:3), each=5000)))
ggplot(data=df) + geom_path( aes(x=step, y=w, col=component),alpha =0.8,size=.8) + theme_bw()+ylim(0,0.8)
dev.off()

