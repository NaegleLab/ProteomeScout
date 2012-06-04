use strict;
use warnings;
use commonTools;
use globalVars;




# handleError($function, $string, $params) 
# Prints error to a log file
# Inputs: $function -name of originating function
#         $string - description of error
#         $params - ref to array of parameter used in function call
# Outputs: writes to log file (appends)
# Kristen Naegle
# March 11, 2008
sub handleError($$$){
    my ($function, $string, $params) = @_;
    my $logFile = returnLogFile();
    open(LOG, ">>$logFile") || die "Can't open log file $logFile for appending\n";
    print LOG "ERROR: $string IN FUNCTION:$function WITH PARAMS: @$params\n";    close(LOG);

}

# handleSysError($function, $string, $params) 
# Prints error to a log file
# Inputs: $function -name of originating function
#         $string - description of error
#         $params - ref to array of parameter used in function call
# Outputs: writes to system log file (appends)
# Kristen Naegle
# September 9, 2009
sub handleSysError($$$){
    my ($function, $string, $params) = @_;
    my $logFile = returnSysLogFile();
    open(LOG, ">>$logFile") || die "Can't open log file $logFile for appending\n";
    my $date = returnTimeStr;
    print LOG "$date\tERROR: $string IN FUNCTION:$function WITH PARAMS: @$params\n";    close(LOG);

}

# handleWarning($function, $string, $params) 
# Prints warning to a log file
# Inputs: $function -name of originating function
#         $string - description of error
#         $params - ref to array of parameter used in function call
# Outputs: writes to log file (appends)
# Kristen Naegle
# March 11, 2008
sub handleWarning($$$){
    my ($function, $string, $params) = @_;
    my $logFile = returnLogFile();
    open(LOG, ">>$logFile") || die "Can't open log file $logFile for appending\n";
    print LOG "WARNING: $string IN FUNCTION:$function WITH PARAMS: @$params\n";    close(LOG);

}

# flushLog($logFile)
# Overwrites to clear log file
# Inputs: $logfile - input file to clear
# Kristen Naegle
# March 11, 2008
sub flushLog(){
    #my $logFile = shift;
    my $logFile = returnLogFile();
    open(LOG, ">$logFile");
    close(LOG);
}

# $logFile = returnLogFile()
# returns location of log file
# Outputs: $logFile - location of logfile
# Kristen Naegle
# March 11, 2008
sub returnLogFile(){
    my $logPath = returnGlobalLogPath();
    
    my $logFile = $logPath."/errorLog";
    if(!-e $logFile) { system("touch $logFile")};
    
    return $logFile;
}

# $processLog = returnProcessLogFile();
# Returns a global processLog file
# Outputs: $processLog - a log file for recording the process name\tprocess params\tstart time\t end time\n.  If end time doesn't exist, check process id. If status killed then know a job was fatally ended.
# Kristen Naegle
# Nov. 25, 2009
sub returnProcessLogFile(){
    my $logPath = returnGlobalLogPath();
    
    my $logFile = $logPath."/processLog";
    if(!-e $logFile) { system("touch $logFile")};
    
    return $logFile;

}

# writeProces($process, $params)
# Writes a line to the processLog with the process, params used, and start time # Inputs: $process - name of process
#         $params - params used in process
# Kristen Naegle
# Nov. 25, 2009
sub writeProcess($$){
    my ($process, $params) = @_;
    my $logFile = returnProcessLogFile();
    open(LOG, ">>$logFile") || die "Can't open log file $logFile for appending\n";

#    print LOG "ERROR: $string IN FUNCTION:$function WITH PARAMS: @$params\n";    close(LOG);
    my $startTime = returnTimeStr();
    print LOG "$process\t@$params\t$startTime\n";
    close(LOG);

}

# writeProcessFinish()
# Append the end time to the last line of the file
# Kristen Naegle
# Nov. 25, 2009
sub writeProcessFinish(){
    my $endTime = returnTimeStr();
    my $logFile = returnProcessLogFile();
    open(LOG, $logFile) || die "Can't open log file $logFile for reading\n";
    my @file = <LOG>;
    close(LOG);
    
    open(LOG, ">$logFile") || die "Can't open log file $logFile for writing\n";
    my $lastLine = $file[$#file-1];
    chomp $lastLine;
    $lastLine .= "\t".$endTime;
    #print "DEBUG last line adding: $lastLine\n";
    $file[$#file-1] = $lastLine;
    $file[$#file] = '';
    print LOG @file;
    close(LOG);
    

}

# $KILLED = checkForKilledProcess
# Checks the process log for a missing end time on the last process
# Ouputs: $KILLED - returns 1 if no end time exists for last process
# Kristen Naegle
# Nov. 25, 2009
sub checkForKilledProcess(){
    my $logFile = returnProcessLogFile();
    open(LOG, $logFile) || die "Can't open log file $logFile for writing\n";
    my @file = <LOG>;
    my $lastLine = $file[$#file-1];
    my @line = split('\t', $lastLine);
    my $KILLED = 0;
    if(scalar(@line) < 4){
	$KILLED = 1;

    }
   # print "DEBUG: line is @line has scalar".scalar(@line)."\n";
    close(LOG);
    return $KILLED;
}

# $logFile = returnSysLogFile()
# returns location of log file
# Outputs: $logFile - location of logfile
# Kristen Naegle
# Sept 10, 2009
sub returnSysLogFile(){
    my $logPath = returnGlobalLogPath();
    my $logFile = $logPath."/system/errorLog";
    if(!-e $logFile) { system("touch $logFile")};
    
    return $logFile;
}

# $logPath = returnGlobalLogPath();
# Returns the global log path
# Outputs: $path - path to (dir)/log
# Kristen Naegle
# Sept 10, 2009
sub returnGlobalLogPath(){
#    my $path = "/home/knaegle/SVN/knaegle/scripts/log";
    my $path = $globalVars::LOG_PATH;
    return $path;
}

# $newLog = dateStampLogFile()
# Moves logfile to a date stamped version for keeping
# Outputs: $newLog - location of date_stamped log file (logfile_year_month_day_hours_min_seconds
# Kristen Naegle
# March 11, 2008
sub dateStampLogFile(){
    my $logFile = returnLogFile();
    my $dateStr = returnTimeStr();
    $dateStr =~ s/:/_/g;
    $dateStr =~ s/-/_/g;
    $dateStr =~ s/ /_/g;
    my $newLog = $logFile."_".$dateStr;
    `cp $logFile $newLog`;
    return $newLog;

}

# printLineToLogFile($);
# Prints a line to the log file
# Inputs: $line - line to print (appends new line feed if doesn't exist)
# Kristen Naegle
# March 18, 2008
sub printLineToLog($){
    my ($line) = @_;
    my $logFile = returnLogFile();
    open(LOG, ">>$logFile") || die "Can't open log file $logFile for appending\n";
    my $temp = $line;
    if ($line eq chomp($temp)){
	$line = $line."\n";
    }
    print LOG $line;
    close(LOG);

}

# printBreakerLineToLog();
# Prints a series of dashes and new line feed to log file
# Kristen Naegle
# March 18, 2008
sub printBreakerLineToLog(){
    my $line = "--------------------------------------------------------------\n";
    printLineToLog($line);

}


1;
