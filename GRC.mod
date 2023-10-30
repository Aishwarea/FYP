NEURON {
	POINT_PROCESS GRC
	RANGE C,gl,el,delT,Vt,I,Vr,a,b,tw,fflag,thresh,gsyn,nspike:,thresh,lastspike
	RANGE tauvt, incr, Vtmax, Isyn, vv,tau
	RANGE i_membrane
}

UNITS {
	(nA)=(nanoamp)
	(pA)=(picoamp)
	(mA)=(microamp)
	(mV)=(millivolt)
	(mM)=(milli/liter)
	(umho)=(micromho)
	(pF)=(picofarad)
	(nS)=(nanosiemens)
}


PARAMETER {
	C=150	(pF) :200
	gl=10	(nS) :10
	el=-70	(mV) :-70
	delT=4	(mV)
	thresh=30	(mV) :0
	Vr=-70	(mV) :-50
	Vt=-50	(mV)
	Isyn=0	(mA)
	a=9	(nS) :2
	b=250	(pA) :0
	tw=20	(ms) :30
	fflag=1
	i_membrane = 0
	tau=9	(ms)
}

ASSIGNED {
	gsyn
	nspike
	lastspike
	t0
}

STATE { 
	w	(pA)
	vv 	(mV)
	:i	(nA)
}

INITIAL {
	vv=-70
	t0=0
	w=0
	i_membrane = 0
	gsyn = 0 :40
	nspike=0
	lastspike=0
	net_send(0,1)
}

BREAKPOINT {
	SOLVE states METHOD derivimplicit

}

DERIVATIVE states {
	Isyn=gsyn*(vv-0)
	vv'=(gl*(el-vv)+gl*delT*exp((vv-Vt)/delT)-Isyn-w)/C
	i_membrane = gl*(el-vv)+gl*delT*exp((vv-Vt)/delT)-Isyn-w
	:vv' = i_membrane/C	
	:printf("vv value is %g\n", vv)
	w'=(a*(vv-el)-w)/tw
	:Vt'=(Vtr-Vt)/tvt
	:gcal'=-gcal/taucal
}

NET_RECEIVE (u) {
	
	INITIAL {
		nspike=0
	}
	
	if (flag == 1) {
		WATCH (vv>thresh) 2
	} else if (flag == 2) {
		net_event(t)
		vv=Vr
		w=w+b
		gsyn=0
		
	} else {	:synaptic activation
		:gsyn = gsyn*exp(-(t - t0)/tau)
		gsyn=gsyn+u
		:t0=t
		printf("gsyn value is %g\n", gsyn)
		printf("u value is %g\n",u)
		
	
	}
}
