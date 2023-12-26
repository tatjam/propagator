#include "Propagator.h"
#include "Output.h"

int main(void)
{
	Propagator prop;
	KeplerElements kepler;
	kepler.a = 7000e3;
	kepler.e = 0.01;
	kepler.raan = 0.4;
	kepler.arg_per = 0.2;
	kepler.inc = 0.4;
	kepler.true_anom = 0.3;
	EulerElements<true> start = kepler_to_euler<true>(kepler);
	prop.init(0.0, start);

	auto results = prop.propagate<true, true>(10000.0, 1.0, 10.0);

	clear_file("out.txt");
	append_table(results, "out.txt");


}