

template<class TA, class TB, class T>
std::pair<result_type<T> > EvaluateWithGuardBits(
`	const TA &a,
	const TB &b,
	int resultWidth
){
	if(a.size() < b.size())
		return MaxErrorWithGuardBits(b,a,f);
	
	result_type<T> resultHU(resultWidth);
	result_type<T> resultHD(resultWidth);
	
	typename TA::const_iterator curr_a=a.begin() end_a=a.end();
	while(curr_a!=end_a){
		T exact_a=curr_a->first, approx_a=curr_a->second;
		
		typename TB::const_iterator curr_b=b.begin() end_b=b.end();
		while(curr_b!=end_b){
			T exact_a=curr_b->first, approx_b=curr_b->second;
			
			T exact=exact_a*exact_b;
			T approx=approx_a*approx_b;
			resultHD.Add(RoundHalfDown(approx));
			resultHU.Add(RoundHalfUp(approx));
			
			curr_b++;
		}
		curr_a++;
	}
}
