 COMMENT
 Last modified on 29-MAR-2014
 Written by: Krishna Chaitanya

 Ampa dynamics taken from Chapter 10 of NEURON book with little bit changes with respect to nspike
 For NMDA dynamics, Mg block is added which is adapted from Migliore et al., 2010 (ModelDB accession: 127995)
 After hyperpolarization dynamics were added to make AdEx model fire properly (adapted from Clopath et al., 2010) 

 References:
 Claudia Clopath and Wulfram Gerstner: Voltage and spike timing interact in STDP – a unified model, Frontiers in Synaptic Neuroscience, 2010. 
 Migliore M, Hines M, McTavish TS, Shepherd GM: Functional roles of distributed synaptic clusters in the mitral-granule cell network of the olfactory  	bulb Frontiers in Integrative Neuroscience 4:122, 2010. 
 ENDCOMMENT

NEURON {
	POINT_PROCESS PURKINJEADEX
	RANGE R, g, gmax, mg, inmda, iampa, igaba, gnmda, gampa
        NONSPECIFIC_CURRENT i
        RANGE Cdur, Alpha, Beta, Rinf, Rtau, nmdafactor
	RANGE r1FIX,r2,r3,r4,r5,r1,r6,r6FIX,kB
	RANGE Ca,gl,el,delT,Vt,I,Vr,a,b,tw,taucal,tvt,tauvr,fflag,thresh,gsyn,nspike:,thresh
	RANGE tauvt, incr, Vtmax, Isyn, v0_block,k_block
	RANGE tau1, tau2, e
}

UNITS {
	(nA)=(nanoamp)
	(pA)=(picoamp)
	(mV)=(millivolt)
	(mM)=(milli/liter)
	(umho)=(micromho)
	(pF)=(picofarad)
	(nS)=(nanosiemens)
	(uS) = (microsiemens)
}


PARAMETER {
	r1FIX		= 5.4		(/ms/mM) 	 
	r2		= 0.82		(/ms)		 
	r3		= 0		(/ms)		 
	r4		= 0		(/ms)		 
	r5		= 0.013		(/ms)		
	r6FIX		= 1.12		(/ms/mM)
	kB		= 0.44		(mM)

	nmdafactor = 0.0035 (1)
	mg	= 1    (mM)		: external magnesium concentration
	gmax = 50200 (umho)
	v0_block 	= -20 		(mV)	: -8.69 (mV)	: -18.69 (mV) : -32.7 (mV)
	k_block 	= 13		(mV)

	Cdur=0.3       (ms)
        Alpha=0.94      (/ms)
        Beta=0.18       (/ms)

	tau1=.1 (ms) <1e-9,1e9>
	tau2 = 100 (ms) <1e-9,1e9>
	e=-70	(mV)

	Ca=100	(pF) :200 350
	gl=10	(nS) :12
	el=-65	(mV) :-58
	delT=2	(mV) :7
	thresh=30	(mV) :0
	Vr=-58	(mV) :-50 -64
	:Vt=-50	(mV) :-60
	Isyn=100	(pA)
	a=-13	(nS) :-12 9
	b=260	(pA) :1460 250
	tw=1	(ms) :30 7
	taucal=40	(ms)
	tvt=50	(ms)
	tauvr=-50.4	(mV)
	nspike=99
	fflag=1
}

ASSIGNED {
	r1		(/ms)
	r6		(/ms)
	Trelease	(mM)
	gampa	(umho)
	gsyn
	i 	(nA)
	Rinf
	Rtau	(ms)
	synon
	inmda 		(nA)		: current = gnmda*(v - E)
	iampa 		(nA)		: current = gampa*(v - E)
	igaba		(nA)
	gnmda 		(umho)
	g (uS)
	factor
}

STATE { 
	C
	O
	D
    	Ron
   	Roff
	w	(pA)
	vv 	(mV)
	gcal	(pA)
	Vt	(mV)
	A (uS)
	B (uS)
}

INITIAL {
	LOCAL tp
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor

	C		=	1
	O		=	0
	D		=	0
	Trelease	=	0 	(mM)
	Rinf=Alpha/(Alpha+Beta)
        Rtau=1/(Alpha+Beta)
        synon=0
	vv=-60	(mV)
	gcal=0
	Vt=-50	(mV)
	w=0	(pA)
	gsyn=0
	net_send(0,1)
}

BREAKPOINT {
	SOLVE states METHOD derivimplicit
	gnmda = mgblock(vv)*(Ron + Roff)*gmax*nmdafactor
	inmda = gnmda*(vv - 0)
	g = B - A
	igaba = g*(vv-e)
	gampa=(Ron+Roff)*1(umho)
        iampa=gampa*(vv-0)
	i=iampa+igaba
}



DERIVATIVE states {
	:Isyn=gsyn*(vv-0)
	Ron'=(synon*Rinf-Ron)/Rtau
        Roff'=-Beta*Roff
	vv'=(gl*(el-vv)+gl*delT*exp((vv-Vt)/delT)-(i+w*(0.001))*(1000)+gcal)/Ca
	w'=(a*(vv-el)-w)/tw
	Vt'=(Vr-Vt)/tvt
	gcal'=-gcal/taucal
	A' = -A/tau1
	B' = -B/tau2
}

FUNCTION mgblock(vv(mV)) {
	TABLE 
	DEPEND mg
	FROM -140 TO 80 WITH 1000

	: from Jahr & Stevens

	mgblock = 1 / (1 + exp(0.062 (/mV) * -vv) * (mg / 3.57 (mM)))
	:mgblock = 1 / ( 1 + exp ( - ( vv - v0_block ) / k_block ) )
}

NET_RECEIVE (u (uS),on,nspike,rth,u_rechan,tth (ms)) {
	:printf("While entereing NET_RECEIVE block, value found is %g\n", flag)
	if (flag == 1) {
		:printf("I entered flag 1\n")
		WATCH (vv>thresh) 2
	} else if (flag == 2) {
		:printf("I entered flag 2\n")
		net_event(t)
		vv=Vr
		w=w+b
		gcal=Isyn
		Vt=tauvr
	} else if (u==0.001) {	:synaptic activation
		:A = A + u*factor :for inhibitory synapse adapted from Exp2syn.mod
		:B = B + u*factor : for inhibitory synapse
		if (flag==0) {
			nspike=nspike+100
			if (!on) {
				rth=rth*exp(-Beta*(t-tth))
               	                tth=t
                                on=1
                        	synon=synon+u
: deprecated
:    	                        state_discontinuity(Ron, Ron+rth)
:                               state_discontinuity(Roff, Roff-rth)
                        	Ron = Ron+rth
                        	Roff = Roff-rth
			}
			:printf("Before nspike value is %g and flag value is %g\n", nspike, flag)
			net_send(Cdur, nspike)
			:printf("nspike value is %g and flag value is %g\n", nspike, flag)
		}
	} else if (u==0.01) {
		:printf("heyeyeyye\n")
		u_rechan=0.00036
		A = A + u_rechan*factor :for inhibitory synapse adapted from Exp2syn.mod
		B = B + u_rechan*factor : for inhibitory synapse
	}
:Added extra flag here to make sure that segmentation error wouldn't happen. 
:Below block undoes the changes happened during the on state of neurotransmitter. 
	if (flag > 2) {
		if (flag==nspike) {
			:printf("Hey i entered flag==nspike\n")
			rth=u*Rinf+(rth-u*Rinf)*exp(-(t-tth)/Rtau)
			tth=t
			synon=synon-u
			Ron=Ron-rth
			Roff=Roff+rth
			on=0
		}
	} 
	:printf("After else block flag value is %g\n", flag)
	COMMENT
  	if (flag != 1) { : prevent segmentation error
		printf("Just above flag==nspike\n")
		printf("flag value here is %g and nspike is %g\n", flag, nspike)
        	if (flag == nspike) {
 			:printf("turn transmitter off\n")
                	:rth=u*Rinf+(rth-u*Rinf)*exp(-(t-tth)/Rtau)
                	:tth=t
                	:synon=synon-u
                	:Ron = Ron-rth
                	:Roff = Roff+rth
                	:on=0
        	}
    	}
	ENDCOMMENT
}
