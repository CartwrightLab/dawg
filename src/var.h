// var.h - Copyright (c) 2004-2009 Reed A. Cartwright (all rights reserved)

#ifndef DAWG_VAR_H
#define DAWG_VAR_H

#include "dawg.h"
#include "tree.h"

class MapSsToVar;

// DawgVar is a variant representing the input variables
class DawgVar
{
public:
	typedef std::vector<DawgVar*> Vec;

	enum Type { tyNone, tyBool, tyNumber, tyString, tyVector, tyTree};

	DawgVar() : m_tyType(tyNone) { }
	virtual ~DawgVar();
	
	// Type Routines
	Type GetType() const { return m_tyType; }
	bool IsType(Type ty) const { return m_tyType == ty; }
	void Unset();

	// Variable Retreival / Setting
	static MapSsToVar& GetMap();
	static DawgVar* GetVar(const std::string &ssKey);
	static void		SetVar(const std::string &ssKey, DawgVar* pVar, int nMode = 0);
	static void		ClearMap();
	
	// General Routines

	// Size returns
	//   0 if tyNone
	//   1 if tyBool, tyNumber, tyString, tyTree
	//   z if tyVector, where z is the length of the vector
	Vec::size_type Size();

	// Number Routines
	explicit DawgVar(double dVar) : m_tyType(tyNumber), m_dData(dVar) { }
	double		GetNumber() const { return m_dData; }
	bool		Get(double &rdVar ) const;
	void		Set(double dVar );
	bool		Get(int &rnVar ) const;
	void		Set(int nVar );
	bool		Get(unsigned int &ruVar ) const;
	void		Set(unsigned int uVar );
	
	// Bool Routines
	explicit DawgVar(bool bVar) : m_tyType(tyBool), m_bData(bVar) { }
	bool		GetBool() const { return m_bData; }
	bool		Get(bool &rbVar ) const;
	void		Set(bool bVar );

 	// String Routines
	explicit DawgVar(const std::string &rssVar) : m_tyType(tyString)
		{ m_pssData = new std::string(rssVar); }
	const std::string&	GetString() const { return *m_pssData; }
	bool				Get(std::string &rssVar ) const;
	void				Set(const std::string &rssVar );
	
	// Vector Routines
	explicit DawgVar(const Vec *p) : m_tyType(tyVector)
		{ m_pvData = p; }
	const Vec*	GetVector() const { return m_pvData; }
	bool		Get(const Vec *&rVec) const;
	void		Set(const Vec* rVec);

	DawgVar& GetAt(Vec::size_type uIndex)
		{ return IsType(tyVector) ? *(*m_pvData)[uIndex] : *this; }
	const DawgVar& GetAt(Vec::size_type uIndex) const
		{ return IsType(tyVector) ? *(*m_pvData)[uIndex] : *this; }
	DawgVar& operator[](Vec::size_type uIndex)
		{ return GetAt(uIndex); }
	const DawgVar& operator[](Vec::size_type uIndex) const
		{ return GetAt(uIndex); }
	
	// Tree Routines
	explicit DawgVar(NewickNode* p) : m_tyType(tyTree), m_ptrData(p) { }
	NewickNode*	GetTree() const { return m_ptrData; }
	bool		Get(NewickNode *&rTree) const;
	void		Set(NewickNode *pTree);

protected:
	Type m_tyType;
	union
	{
		double m_dData;
		bool   m_bData;
		const Vec*   m_pvData;
		const std::string* m_pssData;
		NewickNode*  m_ptrData;
	};
private:
	DawgVar(const DawgVar& var);
	DawgVar& operator = (const DawgVar& var);

public:
	// Templates

	// Get the value of Key and place it in R
	template<class T>
	static bool Get(const std::string &ssKey, T &r)
	{
		DawgVar *pVar = GetVar(ssKey);
		return ( pVar != NULL && pVar->Get(r) );
	}

	// Get as array filling in values as neccessary.
	// An array has a set length.
	template<class T>
	Vec::size_type GetArray(T ar[], Vec::size_type uSize, bool bExpand=true)
	{
		Vec::size_type uMax = std::min(uSize, Size());
		Vec::size_type u = 0;
		// read uMax elements from vector
		for(; u<uMax; u++)
		{
			// stop if element is of wrong type
			if(!GetAt(u).Get(ar[u]))
				return u;
		}
		// fill-in additional elements with the first one
		for(; bExpand && u<uSize; u++)
		{
			ar[u] = ar[0];
		}
		// return number of elements set
		return u;		
	}
	
	// Get as vector.
	// A vector has a variable length.
	template<class T>
	bool GetVector(std::vector<T> &rVec)
	{
		T tTemp;
		rVec.clear();
		// read each element in vector, stoping if of wrong type
		for(Vec::size_type u = 0; u<Size(); u++)
		{
			if(!GetAt(u).Get(tTemp))
				return false;
			rVec.push_back(tTemp);
		}
		return true;
	}

	// Get as matrix, filling in rows as neccessary
	// A matrix is an array of vectors.
	template< class T >
	Vec::size_type GetMatrix(std::vector<T> ar[], Vec::size_type uSize, bool bExpand=true)
	{
		if(Size() == 0)
			return 0;
		Vec::size_type u = 0;
		Vec::size_type uMax = std::min(uSize, Size());

		// Check to see if it is a double vector
		// not 100% accurate
		if(GetAt(0).IsType(tyVector))
		{
			// read rows
			for(;u<uMax;++u)
			{
				if(!GetAt(u).GetVector(ar[u]))
					return u;
			}
		}
		else
		{
			// if single vector, read values into the first row
			if(!GetVector(ar[0]))
				return 0;
			u = 1;
		}
		// fill-in rows as neccessary
		for(; bExpand && u<uSize; u++)
			ar[u] = ar[0];
		return u;
	}
	
	// Get Key as Array
	template< class T >
	static Vec::size_type GetArray( const std::string &ssKey,  T ar[], Vec::size_type uSize, bool bExpand=true)
	{
		DawgVar* pVar = GetVar(ssKey);
		if(pVar == NULL)
			return 0;
		return pVar->GetArray(ar, uSize, bExpand);
	}
	// Get Key as Vector
	template<class T>
	static bool GetVector( const std::string &ssKey, std::vector<T> &rVec)
	{
		DawgVar* pVar = GetVar(ssKey);
		if(pVar == NULL)
			return false;
		return pVar->GetVector(rVec);
	}
	// Get Key as Matrix
	template<class T>
	static Vec::size_type GetMatrix(const std::string &ssKey,  std::vector<T> ar[], Vec::size_type uSize, bool bExpand=true)
	{
		DawgVar* pVar = GetVar(ssKey);
		if(pVar == NULL)
			return 0;
		return pVar->GetMatrix(ar, uSize, bExpand);
	}

};

// A map class that will delete pointers upon destruction
class MapSsToVar : public std::map<std::string, DawgVar*>
{
public:	
	virtual ~MapSsToVar()
	{
		for(iterator it = begin(); it != end(); ++it)
			if(it->second)
				delete it->second;
		clear();
	}
};

#endif //DAWG_VAR_H

