// var.cc - Copyright (C) 2006 Reed A. Cartwright (all rights reserved)

#include "dawg.h"
#include "var.h"

using namespace std;

Variable* VarDB::GetVar(const Key &k)
{
	Map::iterator it = m_map.find(k);
	return (it != m_map.end()) ? it->second : NULL;
}

const Variable* VarDB::GetVar(const Key &k) const
{
	Map::const_iterator it = m_map.find(k);
	return (it != m_map.end()) ? it->second : NULL;
}


// nMode: =, ?=, +=
void VarDB::SetVar(const Key& k, Variable *v, int nMode)
{
	if(nMode == 0)
		m_map[k] = v;
	else if(nMode == 1 && GetVar(k) == NULL)
		m_map[k] = v;
}

// from lexer.ll
bool RunParser(FILE *fin, VarDB *db);

bool VarDB::Parse(const char* cs)
{
	FILE* stream;
	bool bFile = false;
	if(cs == NULL || strcmp(cs, "-") == 0)
	{
		stream = stdin;
		m_csParsedFile = "stdin";
	}
	else
	{
#if _MSC_VER >= 1400
		if(fopen_s(&stream, cs, "r"))
			stream = NULL;
#else	
		stream = fopen(cs, "r");
#endif
		m_csParsedFile = cs;
		bFile = true;
	}
	if(stream == NULL)
		return false;
	bool ret = RunParser(stream, this);
	if(bFile)
		fclose(stream);
	return ret;
}

bool VarDB::ParseError(const char *csMsg, size_t uLine, const char *csText)
{
	cerr << "ALERT: " << csMsg << " in " << m_csParsedFile << " at line " << uLine;
	cerr << ": \"" << csText << "\"." << endl;
	return false;
}
