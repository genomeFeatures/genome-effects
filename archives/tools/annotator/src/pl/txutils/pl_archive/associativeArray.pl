#!/usr/bin/perl

#working with associative arrays
sub dereferenceAss{
 my($array_ref)=@_;
 ${$array_ref}{"a"}="testing";
 ${$array_ref}{"b"}="testingv";
}

%array=("one"=>"1","rwo"=>"1");
dereferenceAss(\%array);
while(($key,$value)=each(%array)){
  print "$key,$value\n";
}

