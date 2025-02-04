#ifndef METHOD_PROGONKY_PQX_H
#define METHOD_PROGONKY_PQX_H

template<typename T>
struct PQ{
	T p_;
	T q_;
	PQ(const T p, const T q) : p_(p), q_(q){}
	PQ() : p_(0), q_(0) {}
};

#endif //METHOD_PROGONKY_PQX_H
