#!/usr/bin/perl
# script to write dependency among files 
# modified by S. Cozzini from the original script:
# Makemake 
# Written by Michael Wester <wester@math.unm.edu> February 16, 1995
# Cotopaxi (Consulting), Albuquerque, New Mexico
#
open(MAKEFILE, "> Makefile.dependencies");
#
#

&MakeDependsf90($ARGV[1]);
#&MakeDepends("*.f",'^\s*include\s+["\']([^"\']+)["\']');
&MakeDepends("*.c",'^\s*#\s*include\s+["\']([^"\']+)["\']');

#
# &PrintWords(current output column, extra tab?, word list); --- print words
#    nicely
#
sub PrintWords {
   local($columns) = 78 - shift(@_);
   local($extratab) = shift(@_);
   local($wordlength);
   #
   print MAKEFILE @_[0];
   $columns -= length(shift(@_));
   foreach $word (@_) {
#      print MAKEFILE " \\\n\t$word";
      $wordlength = length($word);
      if ($wordlength + 1 < $columns) {
         print MAKEFILE " $word";
         $columns -= $wordlength + 1;
         }
      else {
         #
         # Continue onto a new line
         #
         if ($extratab) {
            print MAKEFILE " \\\n\t\t$word";
            $columns = 62 - $wordlength;
            }
         else {
            print MAKEFILE " \\\n\t$word";
            $columns = 70 - $wordlength;
            }
         }
      }
   }


#
# &toLower(string); --- convert string into lower case
#
sub toLower {
   local($string) = @_[0];
   $string =~ tr/A-Z/a-z/;
   $string;
   }

#
# &uniq(sorted word list); --- remove adjacent duplicate words
#
sub uniq {
   local(@words);
   foreach $word (@_) {
      if ($word ne $words[$#words]) {
         push(@words, $word);
         }
      }
   @words;
   }

#
# &MakeDepends(language pattern, include file sed pattern); --- dependency
#    maker
#
sub MakeDepends {
   local(@incs);
   local($lang) = @_[0];
   local($pattern) = @_[1];
#
   foreach $file (<${lang}>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
       print  $file,"\n" ;
      while (<FILE>) { 
         if (/$pattern/i ) {

	     push(@incs, $1) }
         if (/^\s*use\s+([^\s,!]+)/i) {
	     $nom=&toLower($1) ; 
	     $nom="\$\(COMM_LAYER)" if ($nom eq "comm_layer") ;  
	     push(@incs,"\$(OBJDIR)/".$nom.".o") ;
	 }
         }
      if (defined @incs) { 
	  $file =~ s/\.[^.]+$/.o/;
	  print MAKEFILE "\$(OBJDIR)/$file: ";
	  &PrintWords(length($file) + 2, 0,&uniq(sort(@incs)));
	  print MAKEFILE "\n\n";
	  undef @incs;
      }
       undef @incs;
      }
   }

#
# &MakeDependsf90(f90 compiler); --- FORTRAN 90 dependency maker
#
sub MakeDependsf90 {
   local($compiler) = &toLower(@_[0]);
   local(@dependencies);
   local(%filename);
   local(@incs);
   local(@modules);
   local($objfile);
   #
   # Associate each module with the name of the file that contains it
   #
   foreach $file (<*.f90>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /^\s*module\s+([^\s!]+)/i &&
            ($filename{&toLower($1)} = $file) =~ s/\.f90$/.o/;
         }
      }
#
   foreach $file (<*.f>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /^\s*module\s+([^\s!]+)/i &&
            ($filename{&toLower($1)} = $file) =~ s/\.f$/.o/;
         }
      }
   #
   # Print the dependencies of each file that has one or more include's or
   # references one or more modules
   #
   foreach $file (<*.f90>) {
      open(FILE, $file);
      while (<FILE>) {
         /^\s*include\s+["\']([^"\']+)["\']/i && ($nom=$1) ;
          if ( $nom =~ /\.inc$/ ) { push(@incs, $nom); }
         /^\s*use\s+([^\s,!]+)/i && push(@modules, &toLower($1));
         }
      if (defined @incs || defined @modules) {
         ($objfile = $file) =~ s/\.f90$/.o/;
         print MAKEFILE "\$(OBJDIR)/$objfile: ";
         undef @dependencies;
         foreach $module (@modules) {
            push(@dependencies,"\$(OBJDIR)/$filename{$module}");
            }
         @dependencies = &uniq(sort(@dependencies));
         &PrintWords(length($objfile) + 2, 0,
                     @dependencies, &uniq(sort(@incs)));
         print MAKEFILE "\n\n";
         undef @incs;
         undef @modules;
         }
      }
#
   foreach $file (<*.f>) {
      open(FILE, $file);
      while (<FILE>) {
         /^\s*include\s+["\']([^"\']+)["\']/i && ($nom=$1) ;
          if ( $nom =~ /\.inc$/ ) { push(@incs, $nom); }
         /^\s*use\s+([^\s,!]+)/i && push(@modules, &toLower($1));
         }
      if (defined @incs || defined @modules) {
         ($objfile = $file) =~ s/\.f$/.o/;
         print MAKEFILE "\$(OBJDIR)/$objfile: ";
         undef @dependencies;
         foreach $module (@modules) {
            push(@dependencies,"\$(OBJDIR)/$filename{$module}");
            }
         @dependencies = &uniq(sort(@dependencies));
         &PrintWords(length($objfile) + 2, 0,
                     @dependencies, &uniq(sort(@incs)));
         print MAKEFILE "\n\n";
         undef @incs;
         undef @modules;
         }
      }

   }
