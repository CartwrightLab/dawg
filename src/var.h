// var.h - Copyright (C) 2006 Reed A. Cartwright (all rights reserved)
#ifndef DAWG_VAR_H
#define DAWG_VAR_H

#include "dawg.h"
#include "tree.h"

class Variable
{
public:
	typedef std::vector< Variable*> vector_t;

	Variable() { }
	virtual ~Variable() { }

	virtual const char* Type() const { return "Null"; }

	virtual size_t Length() const { return 0; }
	virtual Variable* At(size_t /*pos*/) {return this;}
	virtual const Variable* At(size_t /*pos*/) const { return this;}
	virtual void Append(const Variable * /*v*/) {};

	virtual size_t To(int &/*v*/) const { return 0; }
	virtual size_t To(unsigned int &/*v*/) const { return 0; }
	virtual size_t To(double &/*v*/) const { return 0; }
	virtual size_t To(std::string &/*v*/) const { return 0; }
	virtual size_t To(NewickNode* &/*v*/) const { return 0; }
	virtual size_t To(bool &/*v*/) const { return 0; }

	template<class T, size_t sz>
	size_t To(T (&t)[sz]) const
	{
		size_t i;
		for(i = 0; i < sz; ++i)
			if(!At(i)->To(t[i]))
				return i;
		return i;
	}
	template<class T>
	size_t To(std::vector<T> &t) const
	{
		size_t i;
		size_t sz = Length();
		if(sz == 0)
			return 0;
		t.resize(sz);
		for(i = 0; i < sz; ++i)
			if(!At(i)->To(t[i]))
				return i;
		return i;
	}

protected:

};

class ScalarVar : public Variable
{
public:
	ScalarVar() { }
	virtual size_t length() const { return 1; }
	virtual const char* Type() const { return "Scalar"; }

};

class StringVar : public ScalarVar
{
public:
	StringVar(const std::string& ssVal) : m_ssVal(ssVal) { }
	StringVar(const char *cs) : m_ssVal(cs) { }
	virtual size_t To(std::string &v) const { v = m_ssVal; return 1; }
	virtual const char* Type() const { return "String"; }

protected:
	std::string m_ssVal;	
};

class NumberVar : public ScalarVar
{
public:
	NumberVar(double dVal) : m_dVal(dVal) { }	
	virtual size_t To(double &v) const { v = m_dVal; return 1; }
	virtual size_t To(int &v) const { v = static_cast<int>(m_dVal); return 1; }
	virtual size_t To(unsigned int &v) const { v = static_cast<unsigned int>(m_dVal); return 1; }
	virtual size_t To(bool &v) const { v = (static_cast<int>(m_dVal) != 0); return 1; }
	virtual const char* Type() const { return "Number"; }

protected:
	double m_dVal;
};

class BooleanVar : public ScalarVar
{
public:
	BooleanVar(bool bVal) : m_bVal(bVal) { }
	virtual size_t To(bool &v) { v = m_bVal; return 1; }
	virtual const char* Type() const { return "Boolean"; }

protected:
	bool m_bVal;	
};

class TreeVar : public ScalarVar
{
public:
	TreeVar(NewickNode *p) : m_tree(p) { }
	virtual size_t To(NewickNode *&v) { v = m_tree.get(); return 1; }
	virtual const char* Type() const { return "Tree"; }

protected:
	std::auto_ptr<NewickNode> m_tree;
};

template<class T> void do_delete(T t) { delete t;}

class VectorVar : public Variable
{
public:
	VectorVar(const Variable::vector_t &vVal) : m_vVal(vVal) { }
	VectorVar(const Variable::vector_t::value_type &v) : m_vVal(1, v) { }

	virtual const char* Type() const { return "Vector"; }

	virtual ~VectorVar() {
		std::for_each(m_vVal.begin(), m_vVal.end(), do_delete<Variable*>);
	}

	virtual size_t Length() const { return m_vVal.size(); }
	virtual Variable* At(size_t pos) {return m_vVal.at(pos%Length());}
	virtual const Variable* At(size_t pos) const { return m_vVal.at(pos%Length());}
	virtual void Append(Variable *v)
	{
		m_vVal.push_back(v);
	}

protected:
	Variable::vector_t m_vVal;
};

class VarDB
{
public:
	typedef std::string Key;
	typedef std::map<Key, Variable*> Map;

	const Variable* GetVar(const Key &k) const;
	Variable* GetVar(const Key &k);
	
	void SetVar(const Key &k, Variable *v, int nMode = 0);

	template<class T>
	size_t Get(const Key &k, T &v) const {
		const Variable *p = GetVar(k);
		if(p == NULL)
			return 0;
		return p->To(v);
	}

	bool Parse(const char *cs);

protected:
	Map m_map;

	//Parser Varables
	bool ParseError(const char *csMsg, size_t uLine, const char *csText);
	const char *m_csParsedFile;
	friend void yyerror (char *s);
};

#endif //DAWG_VAR_H
