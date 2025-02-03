#ifndef METHOD_PROGONKY_PQX_H
#define METHOD_PROGONKY_PQX_H

template<typename T>
struct PQX{
	T p_;
	T q_;
	PQX(const T p, const T q) : p_(p), q_(q){}
	PQX() : p_(0), q_(0) {}
};

#endif //METHOD_PROGONKY_PQX_H
