
 NEURON {
	POINT_PROCESS NONBASKET
	NONSPECIFIC_CURRENT i
	RANGE C,gl,el,delT,I,Vr,a,b,tw,fflag,thresh,count,gsyn,nspike:,thresh
	RANGE Vt
	RANGE del, dur, amp, i, i_i,dur_i  
	ELECTRODE_CURRENT i 
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
	C=630	(pF) :200
	gl=30	(nS) :10
	el=-70	(mV) :-70
	delT=2	(mV)
	thresh=-10	(mV) :0
	Vr=-58	(mV) :-50
	Vt=-50	(mV)
	I=500	(pA)
	a=2	(nS) :2
	b=0	(pA) :0
	tw=30	(ms) :30
	fflag=1
	count=0
        i_i=200
	del (ms) : iclamp parameters
	dur (ms)	<0,1e9>
	amp (pA)
}



ASSIGNED {
	gsyn	(ns)
	nspike
	i	(nA)
}

STATE { 
	w	(pA)
	vv 	(mV)
}
INITIAL {
	vv=-60
	w=0
	gsyn=0
	nspike=0
	net_send(0,1)
}


BREAKPOINT {
	SOLVE states METHOD derivimplicit
	:i=gsyn*(vv-0)*(0.001)
		at_time(del) : iclamp conditions
	at_time(del+dur)

	if (t < del + dur && t >= del) {
		i = amp  
	}else{
		i = 0
	}

}


DERIVATIVE states {
	vv'=(gl*(el-vv)+gl*delT*exp((vv-Vt)/delT)+I+(i-w*(0.001))*(1000))/C
	w'=(a*(vv-el)-w)/tw
}

NET_RECEIVE (u(nS)) {
	INITIAL {
		nspike=0
	}
	if (flag == 1) {
		WATCH (vv>thresh) 2
	} else if (flag == 2) {
		net_event(t)
		vv=Vr
		w=w+b
		count=count+1
	} else {	:synaptic activation
		gsyn=gsyn+u
	}
}
