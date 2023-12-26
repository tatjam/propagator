#include "Propagator.h"

int main(void)
{
	Propagator prop;
	KeplerElements kepler;
	EulerElements<true> start = kepler_to_euler<true>(kepler);
	prop.init(0.0, start);

	auto results = prop.propagate<true, true>(1000.0, 1.0, 10.0);


}