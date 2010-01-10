# Biogrep: a multi–threaded pattern matcher for large pattern sets

This tool is designed to quickly 
match large sets of patterns against biosequence databases 
and is optimized for multi–processor computers. Biogrep 
uses standard POSIX extended regular expressions and can 
divide the pattern–matching task between a user–specified 
number of processors.



Biogrep is written in the C programming language using the 
GNU regular expression (Hargreaves & Berry, 1992) and 
POSIX threads (pthreads) (Mueller, 1993) libraries. The pro- 
gram reads query patterns from either a plain text file, one– 
per–line, or from a Teiresias–formated pattern file (Rigout- 
sos & Floratos, 1998). These patterns are treated as POSIX 
extended regular expressions and are searched against a user 
supplied file, which can be either a FastA (Pearson & Lip- 
man, 1988) formatted biosequence database or any text file.