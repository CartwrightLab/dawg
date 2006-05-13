require 'yaml'

skel = <<'EOL'
#ifndef DAWG_DAWGVAR_H
#define DAWG_DAWGVAR_H

#include <string>
#include <vector>
#include <fstream>
#include "var.h"

namespace Dawg {

class Variables
{
public:
#{indent(vardecl,1)}
	bool ParseFile(const char* cs)
	{
		ToDefault();
		VarDB db;
		if(!db.Parse(cs))
			return false;
#{indent(varload,2)}		
		return true;
	}
	void ToDefault()
	{
#{indent(vardef,2)}
	}
	bool LoadStringOrFile(std::string &ss)
	{
		if(ss.empty())
			return false;
		std::ifstream iFile(ss.c_str());
		if(!iFile.is_open())
			return false;
		std::getline(iFile, ss, '\\0');
		return true;
	}
	
};

} //namespace Dawg
#endif

EOL

class DVar
	attr_reader :name
	def initialize(name, hash)
		@name = name
		@dims = Array(hash["Dims"]).reject {|e| (0..1) === e}
		@aliases = Array(hash["Alias"])
		@replaces = Array(hash["Replaces"])
		@text = hash["Text"]
		@default = hash["Default"]
	end
	def prefix
		"var"
	end
	def type
		"Void"
	end
	def btype
		"void"
	end
	def etype
		ty = type
		@dims.reverse_each { |d| ty = "Vector of #{ty}s" unless((0..1) === d) }
		ty		
	end
	def cname
		prefix+@name.gsub(/\./, "")
	end
	def ctype
		f = btype()
		e = ''
		@dims.reverse_each { |d|
			if(d < 0)
				f = "std::vector< #{f} >"
			elsif(d > 1)
				e = "[#{d}]#{e}"
			end
		}
		[f,e]
	end
	def cdecl
		t = ctype;
		ret = "#{t[0]} #{cname}#{t[1]}"
	end
	def carray2(i, v)
		ret = ""
		if(v.kind_of?(Array))
			j = -1
			ret = (v.collect{ |e| j+=1; carray2(j, e) }).flatten
			ret = ret.collect{ |e| "[#{i}]#{e}" }
		else
			ret = "[#{i}] = #{v};\n"
		end
		ret
	end
	def carray(v)
		if(v.kind_of?(Array))
			return '{' + (v.collect {|e| carray(e)}).join(', ') + '}'
		end
		v.to_s
	end
	def cload
		names = [@name] | @aliases #| @replaces
		names = names.compact
		names = names.collect { |e| %Q!db.Get("#{e}", #{cname})! }			
		names = names.join("\n\t|| ")
		ret = ""
		if(@default == nil)
			ret << "if(!(#{names}))\n"
			ret << "\treturn DawgError(\"#{name} is not set properly.  It needs to be a '#{etype}'.\");\n"
		else
			ret << "(#{names});\n"
		end
		ret
	end
	def cprocess(s)
		return s unless(s.kind_of?(String))
		h = {"T" => ctype()[0], "A" => ctype()[1], "$" => cname(), "X" => "" }
		s.gsub(/\$([TAX$])/) {h[$1] }
	end
	def cdefault
		if(@default == nil)
			return ''
		elsif(@default.kind_of?(Array))
			#j = -1
			#v = (@default.collect {|e| j+=1; carray(j, e) }).flatten
			#v = v.collect { |e| "#{cname}#{e}" }
			#return v.join('')
			t = ctype;
			return "static #{t[0]} #{cname}Default#{t[1]} = #{carray(@default)};\nmemcpy(#{cname}, #{cname}Default, sizeof(#{cname}));\n"
		elsif(@default =~ /^\$X/)
			return cprocess(@default)
		else
			return "#{cname} = #{cprocess(@default)};\n"		
		end
		''
	end
	
end

class DBoolean<DVar
	def prefix
		"b"
	end
	def btype
		"bool"
	end
	def type
		"Boolean"
	end
end

class DTree<DVar
	def prefix
		"p"
	end
	def btype
		"NewickNode *"
	end
	def type
		"Tree"
	end
end

class DNumber<DVar
	def prefix
		"d"
	end
	def btype
		"double"
	end
	def type
		"Number"
	end
end

class DNumberInt<DNumber
	def prefix
		"n"
	end
	def btype
		"int"
	end
	def type
		"whole Number"
	end	
end

class DNumberUint<DNumber
	def prefix
		"u"
	end
	def btype
		"unsigned int"
	end
	def type
		"positive whole Number"
	end		
end

class DString<DVar
	def prefix
		"ss"
	end
	def btype
		"std::string"
	end
	def type
		"String"
	end
	def cdefault
		if(@default == nil)
			return ''
		elsif(@default.kind_of?(Array))
			j = -1
			v = (@default.collect {|e| j+=1; carray2(j, e) }).flatten
			v = v.collect { |e| "#{cname}#{e}" }
			return v.join('')
		elsif(@default =~ /^\$X/)
			return cprocess(@default)
		else
			return "#{cname} = #{cprocess(@default)};\n"		
		end
		''
	end	
end

class DStringOrFile < DString
	def cload
		cl = super
		cl << "LoadStringOrFile(#{cname});\n"
	end
	
end

def indent(s, n)
	s.gsub(/(.+(\n|$)|\n)/) { "\t" * n + $1}
end

data = YAML::load($stdin)

$, = "\t"
$\ = "\n"

type = {
	"Number" => proc{|e, h| DNumber.new(e, h)},
	"NumberInt" => proc{|e, h| DNumberInt.new(e, h)},
	"NumberUint" => proc{|e, h| DNumberUint.new(e, h)},
	"String" => proc{|e, h| DString.new(e, h)},
	"StringOrFile" => proc{|e, h| DStringOrFile.new(e, h)},
	"Boolean" => proc{|e, h| DBoolean.new(e, h)},
	"Tree" => proc{|e, h| DTree.new(e, h)},
}

vardecl = ''
varpreload = ''
varload = ''
vardef = ''

(data.sort {|a,b| a[1]['Weight'].to_f <=> b[1]['Weight'].to_f}).each {| e |
	#puts e[0]
	temp = type[e[1]["Type"]].call(e[0], e[1])
	vardecl << temp.cdecl << ";\n"
	varload << temp.cload
	vardef  << temp.cdefault
}



print eval( '%Q(' + skel + ')')
