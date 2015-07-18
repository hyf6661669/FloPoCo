#ifndef __LEXER_CONTEXT_HPP__
#define __LEXER_CONTEXT_HPP__

#include <iostream>
#include <vector>

//using namespace std;

class LexerContext {
public:
	void* scanner;
	int result;
	std::istream* is;
	std::ostream* os;
	int yyTheCycle;
	std::vector<std::pair<std::string, int> > theUseTable;

public:
	LexerContext(std::istream* is = &std::cin, std::ostream* os = &std::cout) {
		init_scanner();
		this->is = is;
		this->os = os;
		yyTheCycle=0;
	}

	//these methods are generated in VHDLLexer.cpp 

	void lex();

	virtual ~LexerContext() { destroy_scanner();}

protected:
	void init_scanner();

	void destroy_scanner();
};



#endif
