use strict;
use warnings;
use DBI;
use globalVars;

# $dbh = returnDBHNOCommit($DB)
# Takes $DB, a single letter code, and returns the proper database handle with no commit
# Inputs: $DB: 'P'- production, 'T'-test, 'U'-update (and modify for more)
# Outputs: $dbh - database handle
# Kristen Naegle
# Dec 17, 2009
sub returnDBHNOCommit($){
    my($DB) = @_;
    my $dbh;
    if($DB eq 'P'){
	$dbh = returnProductionDBHNOCommit();
    }
    elsif($DB eq 'U'){
	$dbh = returnUpdateDBHNOCommit();

    }
    elsif($DB eq 'T'){
	$dbh = returnTestDBHNOCommit();
    }
    else{
	$dbh = -1;
	exit;
    }
    return $dbh;
}

# $dbh = returnTestDBH()
# Outputs: $dbh - DBI database handler
# Returns a database handler to test, with commit on 
# Kristen Naegle
# March 6, 2008
sub returnTestDBH(){ 
# function - open a db handler
    my $driver = 'mysql';
    my $db = "test";
    my $dsn = "dbi:$driver:database=$db";
    my $user = "";
    my $password = "";
#my $password;
    my $dbh = DBI->connect($dsn, $user, $password, {RaiseError => 1, AutoCommit => 1}); # || die $DBI:errstr; #if AutoCommit => 0, must commit changes to database. 
#my $dbh = DBI->connect($dsn, $user, {RaiseError => 1, AutoCommit => 0}); # || die $DBI:errstr;
    if(not defined($dbh)){
	print "ERROR on db connect\n";
	exit;
    }
 return $dbh;   
}

# $dbh = returnTestDBHNOCommit()
# Outputs: $dbh - DBI database handler
# Returns a database handler to test, with commit OFF (must use $dbh->commit() or $dbh->rollback() 
# Kristen Naegle
# March 6, 2008
sub returnTestDBHNOCommit(){ 
# function - open a db handler
    my $driver = 'mysql';
    my $db = $globalVars::TEST_DB;
    my $dsn = "dbi:$driver:database=$db";
    my $user = $globalVars::TESTDB_USER;
    my $password = $globalVars::TESTDB_PWD;
	
#my $password;
    my $dbh = DBI->connect($dsn, $user, $password, {RaiseError => 1, AutoCommit => 0}); # || die $DBI:errstr; #if AutoCommit => 0, must commit changes to database. 
#my $dbh = DBI->connect($dsn, $user, {RaiseError => 1, AutoCommit => 0}); # || die $DBI:errstr;
    if(not defined($dbh)){
	print "ERROR on db connect\n";
	exit;
    }
 return $dbh;   
}

# $dbh = returnProductionDBHNOCommit()
# Outputs: $dbh - DBI database handler
# Returns a database handler to msxplore, with commit OFF (must use $dbh->commit() or $dbh->rollback() 
# Kristen Naegle
# March 19, 2008
sub returnProductionDBHNOCommit(){ 
# function - open a db handler
    my $driver = 'mysql';
    my $db = $globalVars::MYSQL_DB;
    my $dsn = "dbi:$driver:database=$db";
    my $user = $globalVars::MYSQL_USER; 
    my $password = $globalVars::MYSQL_PWD;
    my $dbh = DBI->connect($dsn, $user, $password, {RaiseError => 1, AutoCommit => 0}); # || die $DBI:errstr; #if AutoCommit => 0, must commit changes to database. 
    if(not defined($dbh)){
	print "ERROR on db connect\n";
	exit;
    }
 return $dbh;   
}

# $dbh = returnUpdateDBHNOCommit()
# Outputs: $dbh - DBI database handler
# Returns a database handler to teh update database, with commit OFF (must use $dbh->commit() or $dbh->rollback() 
# Kristen Naegle
# Nov. 24, 2009
sub returnUpdateDBHNOCommit(){ 
# function - open a db handler
    my $driver = 'mysql';
    my $db = $globalVars::MYSQL_UPDATE;
    my $dsn = "dbi:$driver:database=$db";
    my $user = $globalVars::MYSQL_USER; 
    my $password = $globalVars::MYSQL_PWD;
    my $dbh = DBI->connect($dsn, $user, $password, {RaiseError => 1, AutoCommit => 0}); # || die $DBI:errstr; #if AutoCommit => 0, must commit changes to database. 

    if(not defined($dbh)){
	print "ERROR on db connect\n";
	exit;
    }
 return $dbh;   
}


# $dbh = returnProductionDBHWithCommit()
# Outputs: $dbh - DBI database handler
# Returns a database handler to msxplore, with commit ON
# Kristen Naegle
# April 6, 2008
sub returnProductionDBHWithCommit(){ 
# function - open a db handler
   
    my $driver = 'mysql';
    my $db = $globalVars::MYSQL_DB;
    my $dsn = "dbi:$driver:database=$db";
    my $user = $globalVars::MYSQL_USER; 
    my $password = $globalVars::MYSQL_PWD;

    my $dbh = DBI->connect($dsn, $user, $password, {RaiseError => 1, AutoCommit => 1}); # || die $DBI:errstr; #if AutoCommit => 0, must commit changes to database. 
#my $dbh = DBI->connect($dsn, $user, {RaiseError => 1, AutoCommit => 0}); # || die $DBI:errstr;
    if(not defined($dbh)){
	print "ERROR on db connect\n";
	exit;
    }
 return $dbh;   
}



# $dbh = returnRefSeqDBHNOCommit()
# Outputs: $dbh - DBI database handler
# Returns a database handler to refseq, with commit OFF (must use $dbh->commit() or $dbh->rollback() 
# Kristen Naegle
# April 26, 2008
sub returnRefSeqDBHNOCommit(){ 
# function - open a db handler
    my $driver = 'mysql';
    my $db = $globalVars::REFSEQ_DB;
    my $dsn = "dbi:$driver:database=$db";
    my $user = $globalVars::MYSQL_USER;
    my $password = $globalVars::MYSQL_PWD;
    my $dbh = DBI->connect($dsn, $user, $password, {RaiseError => 1, AutoCommit => 0}); # || die $DBI:errstr; #if AutoCommit => 0, must commit changes to database. 
    if(not defined($dbh)){
	print "ERROR on db connect\n";
	exit;
    }
 return $dbh;   
}


1;
