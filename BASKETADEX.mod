
 NEURON {
	POINT_PROCESS BASKET
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
	C=6000	(pF) :100
	gl=12	(nS) :10
	el=-70	(mV) :-65
	delT=2	(mV)
	thresh=35	(mV) :35
	Vr=-58	(mV) :-47
	Vt=-50	(mV)
	I=350	(pA) 
	a=-8	(nS) :-10
	b=0	(pA) :30
	tw=160	(ms) :90
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
