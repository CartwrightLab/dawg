dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/ax_cflags_gcc_option.html
dnl
AC_DEFUN([AX_CFLAGS_GCC_OPTION_OLD], [dnl
AS_VAR_PUSHDEF([FLAGS],[CFLAGS])dnl
AS_VAR_PUSHDEF([VAR],[ac_cv_cflags_gcc_option_$2])dnl
AC_CACHE_CHECK([m4_ifval($1,$1,FLAGS) for gcc m4_ifval($2,$2,-option)],
VAR,[VAR="no, unknown"
 AC_LANG_PUSH(C)
 ac_save_[]FLAGS="$[]FLAGS"
for ac_arg dnl
in "-pedantic  % m4_ifval($2,$2,-option)"  dnl   GCC
   #
do FLAGS="$ac_save_[]FLAGS "`echo $ac_arg | sed -e 's,%%.*,,' -e 's,%,,'`
   AC_TRY_COMPILE([],[return 0;],
   [VAR=`echo $ac_arg | sed -e 's,.*% *,,'` ; break])
done
 FLAGS="$ac_save_[]FLAGS"
 AC_LANG_POP(C)
])
case ".$VAR" in
     .ok|.ok,*) m4_ifvaln($3,$3) ;;
   .|.no|.no,*) m4_ifvaln($4,$4) ;;
   *) m4_ifvaln($3,$3,[
   if echo " $[]m4_ifval($1,$1,FLAGS) " | grep " $VAR " 2>&1 >/dev/null
   then AC_RUN_LOG([: m4_ifval($1,$1,FLAGS) does contain $VAR])
   else AC_RUN_LOG([: m4_ifval($1,$1,FLAGS)="$m4_ifval($1,$1,FLAGS) $VAR"])
                      m4_ifval($1,$1,FLAGS)="$m4_ifval($1,$1,FLAGS) $VAR"
   fi ]) ;;
esac
AS_VAR_POPDEF([VAR])dnl
AS_VAR_POPDEF([FLAGS])dnl
])


dnl the only difference - the LANG selection... and the default FLAGS

AC_DEFUN([AX_CXXFLAGS_GCC_OPTION_OLD], [dnl
AS_VAR_PUSHDEF([FLAGS],[CXXFLAGS])dnl
AS_VAR_PUSHDEF([VAR],[ac_cv_cxxflags_gcc_option_$2])dnl
AC_CACHE_CHECK([m4_ifval($1,$1,FLAGS) for gcc m4_ifval($2,$2,-option)],
VAR,[VAR="no, unknown"
 AC_LANG_PUSH(C++)
 ac_save_[]FLAGS="$[]FLAGS"
for ac_arg dnl
in "-pedantic  % m4_ifval($2,$2,-option)"  dnl   GCC
   #
do FLAGS="$ac_save_[]FLAGS "`echo $ac_arg | sed -e 's,%%.*,,' -e 's,%,,'`
   AC_TRY_COMPILE([],[return 0;],
   [VAR=`echo $ac_arg | sed -e 's,.*% *,,'` ; break])
done
 FLAGS="$ac_save_[]FLAGS"
 AC_LANG_POP(C++)
])
case ".$VAR" in
     .ok|.ok,*) m4_ifvaln($3,$3) ;;
   .|.no|.no,*) m4_ifvaln($4,$4) ;;
   *) m4_ifvaln($3,$3,[
   if echo " $[]m4_ifval($1,$1,FLAGS) " | grep " $VAR " 2>&1 >/dev/null
   then AC_RUN_LOG([: m4_ifval($1,$1,FLAGS) does contain $VAR])
   else AC_RUN_LOG([: m4_ifval($1,$1,FLAGS)="$m4_ifval($1,$1,FLAGS) $VAR"])
                      m4_ifval($1,$1,FLAGS)="$m4_ifval($1,$1,FLAGS) $VAR"
   fi ]) ;;
esac
AS_VAR_POPDEF([VAR])dnl
AS_VAR_POPDEF([FLAGS])dnl
])

dnl -------------------------------------------------------------------------

AC_DEFUN([AX_CFLAGS_GCC_OPTION_NEW], [dnl
AS_VAR_PUSHDEF([FLAGS],[CFLAGS])dnl
AS_VAR_PUSHDEF([VAR],[ac_cv_cflags_gcc_option_$1])dnl
AC_CACHE_CHECK([m4_ifval($2,$2,FLAGS) for gcc m4_ifval($1,$1,-option)],
VAR,[VAR="no, unknown"
 AC_LANG_PUSH(C)
 ac_save_[]FLAGS="$[]FLAGS"
for ac_arg dnl
in "-pedantic  % m4_ifval($1,$1,-option)"  dnl   GCC
   #
do FLAGS="$ac_save_[]FLAGS "`echo $ac_arg | sed -e 's,%%.*,,' -e 's,%,,'`
   AC_TRY_COMPILE([],[return 0;],
   [VAR=`echo $ac_arg | sed -e 's,.*% *,,'` ; break])
done
 FLAGS="$ac_save_[]FLAGS"
 AC_LANG_POP(C)
])
case ".$VAR" in
     .ok|.ok,*) m4_ifvaln($3,$3) ;;
   .|.no|.no,*) m4_ifvaln($4,$4) ;;
   *) m4_ifvaln($3,$3,[
   if echo " $[]m4_ifval($2,$2,FLAGS) " | grep " $VAR " 2>&1 >/dev/null
   then AC_RUN_LOG([: m4_ifval($2,$2,FLAGS) does contain $VAR])
   else AC_RUN_LOG([: m4_ifval($2,$2,FLAGS)="$m4_ifval($2,$2,FLAGS) $VAR"])
                      m4_ifval($2,$2,FLAGS)="$m4_ifval($2,$2,FLAGS) $VAR"
   fi ]) ;;
esac
AS_VAR_POPDEF([VAR])dnl
AS_VAR_POPDEF([FLAGS])dnl
])


dnl the only difference - the LANG selection... and the default FLAGS

AC_DEFUN([AX_CXXFLAGS_GCC_OPTION_NEW], [dnl
AS_VAR_PUSHDEF([FLAGS],[CXXFLAGS])dnl
AS_VAR_PUSHDEF([VAR],[ac_cv_cxxflags_gcc_option_$1])dnl
AC_CACHE_CHECK([m4_ifval($2,$2,FLAGS) for gcc m4_ifval($1,$1,-option)],
VAR,[VAR="no, unknown"
 AC_LANG_PUSH(C++)
 ac_save_[]FLAGS="$[]FLAGS"
for ac_arg dnl
in "-pedantic  % m4_ifval($1,$1,-option)"  dnl   GCC
   #
do FLAGS="$ac_save_[]FLAGS "`echo $ac_arg | sed -e 's,%%.*,,' -e 's,%,,'`
   AC_TRY_COMPILE([],[return 0;],
   [VAR=`echo $ac_arg | sed -e 's,.*% *,,'` ; break])
done
 FLAGS="$ac_save_[]FLAGS"
 AC_LANG_POP(C++)
])
case ".$VAR" in
     .ok|.ok,*) m4_ifvaln($3,$3) ;;
   .|.no|.no,*) m4_ifvaln($4,$4) ;;
   *) m4_ifvaln($3,$3,[
   if echo " $[]m4_ifval($2,$2,FLAGS) " | grep " $VAR " 2>&1 >/dev/null
   then AC_RUN_LOG([: m4_ifval($2,$2,FLAGS) does contain $VAR])
   else AC_RUN_LOG([: m4_ifval($2,$2,FLAGS)="$m4_ifval($2,$2,FLAGS) $VAR"])
                      m4_ifval($2,$2,FLAGS)="$m4_ifval($2,$2,FLAGS) $VAR"
   fi ]) ;;
esac
AS_VAR_POPDEF([VAR])dnl
AS_VAR_POPDEF([FLAGS])dnl
])

AC_DEFUN([AX_CFLAGS_GCC_OPTION],[ifelse(m4_bregexp([$2],[-]),-1,
[AX_CFLAGS_GCC_OPTION_NEW($@)],[AX_CFLAGS_GCC_OPTION_OLD($@)])])

AC_DEFUN([AX_CXXFLAGS_GCC_OPTION],[ifelse(m4_bregexp([$2],[-]),-1,
[AX_CXXFLAGS_GCC_OPTION_NEW($@)],[AX_CXXFLAGS_GCC_OPTION_OLD($@)])])

dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/ac_cxx_member_templates.html
dnl
AC_DEFUN([AC_CXX_MEMBER_TEMPLATES],
[AC_CACHE_CHECK(whether the compiler supports member templates,
ac_cv_cxx_member_templates,
[AC_LANG_PUSH(C++)
 AC_TRY_COMPILE([
template<class T, int N> class A
{ public:
  template<int N2> A<T,N> operator=(const A<T,N2>& z) { return A<T,N>(); }
};],[A<double,4> x; A<double,7> y; x = y; return 0;],
 ac_cv_cxx_member_templates=yes, ac_cv_cxx_member_templates=no)
 AC_LANG_POP(C++)
])
if test "$ac_cv_cxx_member_templates" = yes; then
  AC_DEFINE(HAVE_MEMBER_TEMPLATES,,[define if the compiler supports member templates])
fi
])

dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/ac_cxx_member_templates_outside_class.html
dnl
AC_DEFUN([AC_CXX_MEMBER_TEMPLATES_OUTSIDE_CLASS],
[AC_CACHE_CHECK(whether the compiler supports member templates outside the class declaration,
ac_cv_cxx_member_templates_outside_class,
[AC_LANG_PUSH(C++)
 AC_TRY_COMPILE([
template<class T, int N> class A
{ public :
  template<int N2> A<T,N> operator=(const A<T,N2>& z);
};
template<class T, int N> template<int N2>
A<T,N> A<T,N>::operator=(const A<T,N2>& z){ return A<T,N>(); }],[
A<double,4> x; A<double,7> y; x = y; return 0;],
 ac_cv_cxx_member_templates_outside_class=yes, ac_cv_cxx_member_templates_outside_class=no)
 AC_LANG_POP(C++)
])
if test "$ac_cv_cxx_member_templates_outside_class" = yes; then
  AC_DEFINE(HAVE_MEMBER_TEMPLATES_OUTSIDE_CLASS,,
            [define if the compiler supports member templates outside the class declaration])
fi
])
