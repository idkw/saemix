############################### Simulation - MCMC kernels (E-step) #############################

# Helper: compute log-likelihood with IOV
# For IOV, phi at observation level = mean.phi[subject,] + eta[subject,] + gamma[idocc,iov_cols] + beta.occ[occ,]
# We build phiM.obs from the components and pass to compute.LLy
compute.LLy.iov<-function(phiM, gammaM, beta.occ, Uargs, Dargs, DYF, pres) {
	# Build observation-level phi including IOV effects
	# phiM is NM x nb.parameters (subject level, with eta already included)
	# gammaM is NM.occ x nb.iovas (IOV random effects per subject-occasion pair)
	# beta.occ is nocc x nb.parameters (occasion fixed effects)

	# Expand phiM to observation level via IdM, then add gamma and beta.occ
	nb.parameters<-Uargs$nb.parameters
	nobs.total<-length(Dargs$yM)

	# phi at observation level = phiM[subject_of_obs, ]
	phiM.obs<-phiM[Dargs$IdM,]

	# Add gamma (IOV random effects) for IOV parameters
	i1.iov<-Uargs$i1.iov
	gamma.obs<-matrix(0, nrow=nobs.total, ncol=nb.parameters)
	gamma.obs[, i1.iov]<-gammaM[Dargs$idocc.of.obsM, , drop=FALSE]
	phiM.obs<-phiM.obs + gamma.obs

	# Add beta.occ (occasion fixed effects)
	beta.obs<-beta.occ[Dargs$occM, , drop=FALSE]
	phiM.obs<-phiM.obs + beta.obs

	# Transform to psi and compute predictions
	psiM.obs<-transphi(phiM.obs, Dargs$transform.par)
	fpred<-Dargs$structural.model(psiM.obs, 1:nobs.total, Dargs$XM)
	for(ityp in Dargs$etype.exp) fpred[Dargs$XM$ytype==ityp]<-log(cutoff(fpred[Dargs$XM$ytype==ityp]))

	if (Dargs$modeltype=="structural"){
		gpred<-error(fpred, pres, Dargs$XM$ytype)
		# We need to sum by subject (IdM) to get U per subject chain
		lly.obs<-0.5*((Dargs$yM-fpred)/gpred)**2+log(gpred)
	} else {
		lly.obs<- -fpred
	}
	# Sum by subject (IdM maps each obs to a chain-specific subject index)
	U<-tapply(lly.obs, Dargs$IdM, sum)
	U<-as.numeric(U[as.character(1:Dargs$NM)])
	U[is.na(U)]<-0
	return(U)
}


estep<-function(kiter, Uargs, Dargs, opt, mean.phi, varList, DYF, phiM, gammaM=NULL, beta.occ=NULL) {
	# E-step - simulate unknown parameters
	# Input: kiter, Uargs, mean.phi (unchanged)
	# Output: varList, DYF, phiM (changed), and optionally gammaM for IOV

	# Function to perform MCMC simulation
	nb.etas<-length(varList$ind.eta)
	domega<-cutoff(mydiag(varList$omega[varList$ind.eta,varList$ind.eta]),.Machine$double.eps)
	omega.eta<-varList$omega[varList$ind.eta,varList$ind.eta,drop=FALSE]
	omega.eta<-omega.eta-mydiag(mydiag(varList$omega[varList$ind.eta,varList$ind.eta]))+mydiag(domega)
	chol.omega<-try(chol(omega.eta))
	somega<-solve(omega.eta)

	# "/" dans Matlab = division matricielle, selon la doc "roughly" B*INV(A) (et *= produit matriciel...)

	VK<-rep(c(1:nb.etas),2)
	mean.phiM<-do.call(rbind,rep(list(mean.phi),Uargs$nchains))
	phiM[,varList$ind0.eta]<-mean.phiM[,varList$ind0.eta]

	# IOV setup
	has.iov<-Dargs$has.iov
	if(has.iov) {
		nb.iovas<-Uargs$nb.iovas
		i1.iov<-Uargs$i1.iov
		psi.iov<-varList$psi.iov
		domega.iov<-cutoff(mydiag(psi.iov[i1.iov,i1.iov,drop=FALSE]),.Machine$double.eps)
		psi.iov.eta<-psi.iov[i1.iov,i1.iov,drop=FALSE]
		psi.iov.eta<-psi.iov.eta-mydiag(mydiag(psi.iov.eta))+mydiag(domega.iov)
		chol.psi<-try(chol(psi.iov.eta))
		spsi<-solve(psi.iov.eta)
		gammaMc<-gammaM
	}

	# Compute initial log-likelihood
	if(has.iov) {
		U.y<-compute.LLy.iov(phiM, gammaM, beta.occ, Uargs, Dargs, DYF, varList$pres)
	} else {
		U.y<-compute.LLy(phiM,Uargs,Dargs,DYF,varList$pres)
	}

	etaM<-phiM[,varList$ind.eta]-mean.phiM[,varList$ind.eta,drop=FALSE]
	phiMc<-phiM

	# ---- Kernel 1: Full dimension Gaussian proposal for eta ----
	for(u in 1:opt$nbiter.mcmc[1]) {
		etaMc<-matrix(rnorm(Dargs$NM*nb.etas),ncol=nb.etas)%*%chol.omega
		phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc
		if(has.iov) {
			Uc.y<-compute.LLy.iov(phiMc, gammaM, beta.occ, Uargs, Dargs, DYF, varList$pres)
		} else {
			Uc.y<-compute.LLy(phiMc,Uargs,Dargs,DYF,varList$pres)
		}
		deltau<-Uc.y-U.y
		ind<-which(deltau<(-1)*log(runif(Dargs$NM)))
		etaM[ind,]<-etaMc[ind,]
		U.y[ind]<-Uc.y[ind]
	}
	U.eta<-0.5*rowSums(etaM*(etaM%*%somega))

	# ---- Kernel 2: Univariate random walk for eta ----
	if(opt$nbiter.mcmc[2]>0) {
		nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
		nrs2<-1
		for (u in 1:opt$nbiter.mcmc[2]) {
			for(vk2 in 1:nb.etas) {
				etaMc<-etaM
				etaMc[,vk2]<-etaM[,vk2]+matrix(rnorm(Dargs$NM*nrs2), ncol=nrs2)%*%mydiag(varList$domega2[vk2,nrs2],nrow=1)
				phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc
				if(has.iov) {
					Uc.y<-compute.LLy.iov(phiMc, gammaM, beta.occ, Uargs, Dargs, DYF, varList$pres)
				} else {
					Uc.y<-compute.LLy(phiMc,Uargs,Dargs,DYF,varList$pres)
				}
				Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))
				deltu<-Uc.y-U.y+Uc.eta-U.eta
				ind<-which(deltu<(-1)*log(runif(Dargs$NM)))
				etaM[ind,]<-etaMc[ind,]
				U.y[ind]<-Uc.y[ind]
				U.eta[ind]<-Uc.eta[ind]
				nbc2[vk2]<-nbc2[vk2]+length(ind)
				nt2[vk2]<-nt2[vk2]+Dargs$NM
			}
		}
		varList$domega2[,nrs2]<-varList$domega2[,nrs2]*(1+opt$stepsize.rw* (nbc2/nt2-opt$proba.mcmc))
	}

	# ---- Kernel 3: Grouped dimension updates for eta ----
	if(opt$nbiter.mcmc[3]>0) {
		nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
		nrs2<-kiter%%(nb.etas-1)+2
		if(is.nan(nrs2)) nrs2<-1 # to deal with case nb.etas=1
		for (u in 1:opt$nbiter.mcmc[3]) {
			if(nrs2<nb.etas) {
				vk<-c(0,sample(c(1:(nb.etas-1)),nrs2-1))
				nb.iter2<-nb.etas
			} else {
				vk<-0:(nb.etas-1)
				nb.iter2<-1
			}
			for(k2 in 1:nb.iter2) {
				vk2<-VK[k2+vk]
				etaMc<-etaM
				etaMc[,vk2]<-etaM[,vk2]+matrix(rnorm(Dargs$NM*nrs2), ncol=nrs2)%*%mydiag(varList$domega2[vk2,nrs2])
				phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc
				if(has.iov) {
					Uc.y<-compute.LLy.iov(phiMc, gammaM, beta.occ, Uargs, Dargs, DYF, varList$pres)
				} else {
					Uc.y<-compute.LLy(phiMc,Uargs,Dargs,DYF,varList$pres)
				}
				Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%somega))
				deltu<-Uc.y-U.y+Uc.eta-U.eta
				ind<-which(deltu<(-log(runif(Dargs$NM))))
				etaM[ind,]<-etaMc[ind,]
				U.y[ind]<-Uc.y[ind]
				U.eta[ind]<-Uc.eta[ind]
				nbc2[vk2]<-nbc2[vk2]+length(ind)
				nt2[vk2]<-nt2[vk2]+Dargs$NM
			}
		}
		varList$domega2[,nrs2]<-varList$domega2[,nrs2]*(1+opt$stepsize.rw* (nbc2/nt2-opt$proba.mcmc))
	}

	# Update phiM with current eta
	phiM[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaM

	# ---- IOV Kernel: MCMC for gamma (IOV random effects) ----
	if(has.iov) {
		# Recompute U.y with current phiM
		U.y<-compute.LLy.iov(phiM, gammaM, beta.occ, Uargs, Dargs, DYF, varList$pres)
		# Prior on gamma: U.gamma = 0.5 * sum gamma_ik' Psi^{-1} gamma_ik
		U.gamma<-0.5*rowSums(gammaM*(gammaM%*%spsi))
		# Map U.gamma from (NM.occ) to subjects (NM) by summing
		U.gamma.subj<-tapply(U.gamma, Dargs$id.of.idoccM, sum)
		U.gamma.subj<-as.numeric(U.gamma.subj[as.character(1:Dargs$NM)])
		U.gamma.subj[is.na(U.gamma.subj)]<-0

		# Univariate random walk for each IOV component
		nt2.iov<-nbc2.iov<-matrix(0, nrow=nb.iovas, ncol=1)
		for(u in 1:max(1, opt$nbiter.mcmc[2])) {
			for(vk2 in 1:nb.iovas) {
				gammaMc<-gammaM
				gammaMc[,vk2]<-gammaM[,vk2]+rnorm(Dargs$NM.occ)*varList$domega2.iov[vk2,1]
				Uc.y<-compute.LLy.iov(phiM, gammaMc, beta.occ, Uargs, Dargs, DYF, varList$pres)
				Uc.gamma<-0.5*rowSums(gammaMc*(gammaMc%*%spsi))
				Uc.gamma.subj<-tapply(Uc.gamma, Dargs$id.of.idoccM, sum)
				Uc.gamma.subj<-as.numeric(Uc.gamma.subj[as.character(1:Dargs$NM)])
				Uc.gamma.subj[is.na(Uc.gamma.subj)]<-0
				deltu<-Uc.y-U.y+Uc.gamma.subj-U.gamma.subj
				ind<-which(deltu<(-1)*log(runif(Dargs$NM)))
				# Accept: update gammaM for all (subject,occasion) pairs of accepted subjects
				if(length(ind)>0) {
					accepted.idocc<-which(Dargs$id.of.idoccM %in% ind)
					gammaM[accepted.idocc,]<-gammaMc[accepted.idocc,]
					U.y[ind]<-Uc.y[ind]
					U.gamma.subj[ind]<-Uc.gamma.subj[ind]
				}
				nbc2.iov[vk2]<-nbc2.iov[vk2]+length(ind)
				nt2.iov[vk2]<-nt2.iov[vk2]+Dargs$NM
			}
		}
		varList$domega2.iov[,1]<-varList$domega2.iov[,1]*(1+opt$stepsize.rw*(nbc2.iov/nt2.iov-opt$proba.mcmc))
	}


	# ---- Kernel 4: MAP-based proposal (skip for IOV for now) ----
	if(opt$nbiter.mcmc[4]>0 & kiter<opt$nbiter.map & !has.iov) {
		etaMc<-etaM
		propc <- U.eta
		prop <- U.eta
		phi.map<-mean.phi
	  	i1.omega2<-varList$ind.eta
	    iomega.phi1<-solve(omega.eta[i1.omega2,i1.omega2])

	 	# Setup for MAP calculation (MAP is identical no matter the chain)
	  	id<-Dargs$IdM[1:Dargs$nobs]
	  	xind<-Dargs$XM[1:Dargs$nobs,]
	  	yobs<-Dargs$yM[1:Dargs$nobs]
	  	id.list<-unique(id)

	  	if(Dargs$modeltype=="structural"){
	  		for(i in 1:length(id.list)) {
  	  	  		isuj<-id.list[i]
			    xi<-xind[id==isuj,,drop=FALSE]
			    yi<-yobs[id==isuj]
			    idi<-rep(1,length(yi))
			    mean.phi1<-mean.phiM[i,i1.omega2]
			    phii<-phiM[i,]
			    phi1<-phii[i1.omega2]
			    suppressWarnings(phi1.opti<-optim(par=phi1, fn=conditional.distribution_c, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=Dargs$transform.par, model=Dargs$structural.model, pres=varList$pres, err=Dargs$error.model))
			    phi.map[i,i1.omega2]<-phi1.opti$par
			}

			# Repeat the map nchains time
			phi.map <- phi.map[rep(seq_len(nrow(phi.map)),Uargs$nchains ), ]

		 	map.psi<-transphi(phi.map,Dargs$transform.par)
			map.psi<-data.frame(id=id.list,map.psi)
			map.phi<-data.frame(id=id.list,phi.map)
			psi_map <- as.matrix(map.psi[,-c(1)])
			phi_map <- as.matrix(map.phi[,-c(1)])
			eta_map <- phi_map - mean.phiM

			fpred1<-Dargs$structural.model(psi_map, Dargs$IdM, Dargs$XM)
			gradf <- matrix(0L, nrow = length(fpred1), ncol = nb.etas)

			## Compute gradient of structural model (gradf)
			for (j in 1:nb.etas) {
				psi_map2 <- psi_map
				psi_map2[,j] <- psi_map[,j]+psi_map[,j]/1000
				fpred1<-Dargs$structural.model(psi_map, Dargs$IdM, Dargs$XM)
				fpred2<-Dargs$structural.model(psi_map2, Dargs$IdM, Dargs$XM)
				for (i in 1:(Dargs$NM)){
					r = which(Dargs$IdM==i)
					gradf[r,j] <- (fpred2[r] - fpred1[r])/(psi_map[i,j]/1000)
				}
			}


			## Compute gradient of mapping (psi to phi) function (gradh)
			gradh <- list(omega.eta,omega.eta)
			for (i in 1:Dargs$NM){
				gradh[[i]] <- gradh[[1]]
			}
			for (j in 1:nb.etas) {
				phi_map2 <- phi_map
				phi_map2[,j] <- phi_map[,j]+phi_map[,j]/1000
				psi_map2 <- transphi(phi_map2,Dargs$transform.par)
				for (i in 1:(Dargs$NM)){
					gradh[[i]][,j] <- (psi_map2[i,] - psi_map[i,])/(phi_map[i,]/1000)
				}
			}

			## Calculation of the covariance matrix of the proposal
			Gamma <- chol.Gamma <- inv.chol.Gamma <- inv.Gamma <- list(omega.eta,omega.eta)
			for (i in 1:(Dargs$NM)){
				r = which(Dargs$IdM==i)
		        temp <- gradf[r,]%*%gradh[[i]]
				Gamma[[i]] <- solve(t(temp)%*%temp/(varList$pres[1])^2+solve(omega.eta))
				chol.Gamma[[i]] <- chol(Gamma[[i]])
				inv.chol.Gamma[[i]] <- solve(chol.Gamma[[i]])
				inv.Gamma[[i]] <- solve(Gamma[[i]])
			}

		} else {
		  for(i in 1:length(id.list)) {
			    isuj<-id.list[i]
			    xi<-xind[id==isuj,,drop=FALSE]
			    yi<-yobs[id==isuj]
			    idi<-rep(1,length(yi))
			    mean.phi1<-mean.phiM[i,i1.omega2]
			    phii<-phiM[i,]
			    phi1<-phii[i1.omega2]
			    suppressWarnings(phi1.opti<-optim(par=phi1, fn=conditional.distribution_d, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=Dargs$transform.par, model=Dargs$structural.model))
			    phi.map[i,i1.omega2]<-phi1.opti$par
			}
			#rep the map nchains time
			phi.map <- phi.map[rep(seq_len(nrow(phi.map)),Uargs$nchains ), ]
		  	map.psi<-transphi(phi.map,Dargs$transform.par)
			map.psi<-data.frame(id=id.list,map.psi)
			map.phi<-data.frame(id=id.list,phi.map)

			psi_map <- as.matrix(map.psi[,-c(1)])
			phi_map <- as.matrix(map.phi[,-c(1)])
			eta_map <- phi_map[,varList$ind.eta] - mean.phiM[,varList$ind.eta]

			#gradient at the map estimation
			gradp <- matrix(0L, nrow = Dargs$NM, ncol = nb.etas)

			for (j in 1:nb.etas) {
				phi_map2 <- phi_map
				phi_map2[,j] <- phi_map[,j]+phi_map[,j]/100;
				psi_map2 <- transphi(phi_map2,Dargs$transform.par)
				fpred1<-Dargs$structural.model(psi_map, Dargs$IdM, Dargs$XM)
				DYF[Uargs$ind.ioM]<- fpred1
				l1<-colSums(DYF)
				fpred2<-Dargs$structural.model(psi_map2, Dargs$IdM, Dargs$XM)
				DYF[Uargs$ind.ioM]<- fpred2
				l2<-colSums(DYF)

				for (i in 1:(Dargs$NM)){
					gradp[i,j] <- (l2[i] - l1[i])/(phi_map[i,j]/100)
				}
			}

			#calculation of the covariance matrix of the proposal
			fpred<-Dargs$structural.model(psi_map, Dargs$IdM, Dargs$XM)
			DYF[Uargs$ind.ioM]<- fpred
			denom <- colSums(DYF)

			Gamma <- chol.Gamma <- inv.Gamma <- list(omega.eta,omega.eta)
			z <- matrix(0L, nrow = length(fpred), ncol = 1)
			for (i in 1:(Dargs$NM)){
				Gamma[[i]] <- solve(gradp[i,]%*%t(gradp[i,])/denom[i]^2+solve(omega.eta))
				chol.Gamma[[i]] <- chol(Gamma[[i]])
				inv.Gamma[[i]] <- solve(Gamma[[i]])
			}

		}

		etaM <- eta_map
	  	phiM<-etaM+mean.phiM
	  	U.eta<-0.5*rowSums(etaM*(etaM%*%somega))
	  	U.y<-compute.LLy(phiM,Uargs,Dargs,DYF,varList$pres)

	  	for (u in 1:opt$nbiter.mcmc[4]) {

			#generate candidate eta with new proposal
			for (i in 1:(Dargs$NM)){
				Mi <- rnorm(nb.etas)%*%chol.Gamma[[i]]
				etaMc[i,varList$ind.eta]<- eta_map[i,varList$ind.eta] + Mi
			}

			phiMc[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaMc[,varList$ind.eta]
			Uc.y<-compute.LLy(phiM,Uargs,Dargs,DYF,varList$pres)
			Uc.eta<-0.5*rowSums(etaMc[,varList$ind.eta]*(etaMc[,varList$ind.eta]%*%somega))

			for (i in 1:(Dargs$NM)){
				propc[i] <- 0.5*rowSums((etaMc[i,varList$ind.eta]-eta_map[i,varList$ind.eta])*(etaMc[i,varList$ind.eta]-eta_map[i,varList$ind.eta])%*%inv.Gamma[[i]])
				prop[i] <- 0.5*rowSums((etaM[i,varList$ind.eta]-eta_map[i,varList$ind.eta])*(etaM[i,varList$ind.eta]-eta_map[i,varList$ind.eta])%*%inv.Gamma[[i]])
			}

			deltu<-Uc.y-U.y+Uc.eta-U.eta + prop - propc
			ind<-which(deltu<(-1)*log(runif(Dargs$NM)))
			etaM[ind,varList$ind.eta]<-etaMc[ind,varList$ind.eta]
			U.y[ind]<-Uc.y[ind]
			U.eta[ind]<-Uc.eta[ind]

  		}

		phiM[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaM
	}

	# Final phiM update (for non-MAP path)
	if(!(opt$nbiter.mcmc[4]>0 & kiter<opt$nbiter.map & !has.iov)) {
		phiM[,varList$ind.eta]<-mean.phiM[,varList$ind.eta]+etaM
	}

	return(list(varList=varList,DYF=DYF,phiM=phiM, etaM=etaM, gammaM=gammaM))
}
