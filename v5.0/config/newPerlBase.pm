package newPerlBase;
use Exporter;
use strict;
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
use File::Path qw(make_path);
use FindBin qw($Bin $Script);

BEGIN {
        our @ISA = qw(Exporter);
        our @EXPORT = qw(selectconf readconf logAndDie runOrDie createLog timeLog writeLog totalTime stepTime mkdirOrDie stepStart runOrDieWithThreads);
        our $VERSION = 2.0;
}

our $userlog;
our $syslog;
our $BEGIN_TIME=time();
my $nodename=`hostname`;chomp $nodename;
our %stepstarttime;
our %stepname;
sub createLog{
	my $name=shift;
	my $version=shift;
	my $pid=shift;
	my $userpath=shift;
	my $debug=shift;
	my $user=`whoami`;chomp $user;
	mkdirOrDie($userpath) if(!-d $userpath);
	$userpath=abs_path($userpath);
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime($BEGIN_TIME);
	my $date=sprintf("%4d%02d%02d", $year+1900, $mon+1, $day);
	my $prefix=join("-",$name,$version,$date,$user,$pid,$BEGIN_TIME,$nodename);    #log file name prefix
	$prefix.="-debug" if(defined $debug);
	$userlog="$userpath/$prefix.log";
	system("touch $userlog");
	my $workDir = `pwd`;chomp $workDir;
	my $cmdlog=`ps -p $pid -ocmd --no-header`;chomp $cmdlog;
	my $info="=======================================================\n"."start time:\t".($year+1900)."-".($mon+1)."-$day $hour:$min:$sec\n"."user:\t\t$user\n"."cmd:\t\t$cmdlog\n"."work directory:\t$workDir\n"."=======================================================\n";
	writeLog($info);
}
sub selectconf{
    my $conf=shift;
    my $cfg=readconf($conf);
    return $cfg;
}
sub readconf{
	my $file=shift;
	my %cfg;
	open IN,$file or logAndDie("$file not exists!");
	while (<IN>) {
		chomp;
		next if($_=~/^#/);
		next if($_=~/^\s*$/);
		$_=~s/\r//g;
		$_=~s/^\s+//g;
		$_=~s/\s+$//g;
		my ($key,$value)=split /\s+/,$_,2;
		$cfg{$key}=$value;
	}
	close IN;
	return \%cfg;
}

sub runOrDie{
	my ($cmd) = @_ ;
	if($cmd!~/\s/){
		$cmd=abs_path($cmd);
		$cmd="sh $cmd"; # for shell file
	}
	&timeLog($cmd);
	my $begintime=time();
	my $flag = system($cmd);
	if ($flag != 0){
		&timeLog("Error: command fail: $cmd");
		exit(1);
	}
	infoTime($begintime,"command done!");
}
sub runOrDieWithThreads{
	my ($thread,@cmds) = @_ ;
	my $pm = Parallel::ForkManager->new($thread);
	for my $i (@cmds) {
		my $pid = $pm->start and next;
		runOrDie($i);
		$pm->finish;
	}
	$pm->wait_all_children;
}

sub totalTime {
	logAndDie("Start time not recorded!") if (!defined $BEGIN_TIME);
	infoTime($BEGIN_TIME,"All Analysis Finished!");
}
sub stepStart {
	my $num=shift;
	my $info=shift;
	$stepname{"step$num"}=$info;
	$stepstarttime{"step$num"}=time();
	my $detail="Step_$num: $info: start";
	timeLog($detail);
}
sub infoTime {
	my $starttime=shift;
	my $info=shift;
	my $els=time()-$starttime;
	my $detail="$info Elapsed time: $els s";
	timeLog($detail);
}
sub stepTime {
	my $num=shift;
	my $starttime;
	logAndDie("step_$num not started yet") if (!exists $stepname{"step$num"});
	my $info=$stepname{"step$num"};
	my $els=time()-$stepstarttime{"step$num"};
	my $detail="Step_$num: $info: finished! Elapsed time: $els s";
	timeLog($detail);
}
sub logAndDie {
	my $info = shift;
	my ($package, $filename, $line) = caller;
	if(!defined $info){
		timeLog("logAndDie: info is undefined");
		die;
	}
	timeLog("FATAL ERROR (file $filename line $line): $info");
	die;
}
my $formatDateTime = sub{
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	my $format_time = sprintf("%4d-%02d-%02d %02d:%02d:%02d", 
		$year+1900, $mon+1, $day, $hour, $min, $sec);
	return $format_time;
};

sub timeLog {
	my $detail = shift;
	# get current time with string
	my $curr_time = &$formatDateTime(localtime(time()));
	# print info with time
	writeLog("[$curr_time] $detail");
}
sub writeLog{
	my $detail=shift;
	if(defined $userlog and defined $syslog){
		open OUT1,">>$userlog" or die("$userlog PATH ERROR!\n");
		open OUT2,">>$syslog" or die("$syslog PATH ERROR!\n");
		print OUT1 "$detail \n";
		print OUT2 "$detail \n";
		close OUT1;
		close OUT2;
	}elsif(defined $userlog and !defined $syslog){
		open OUT1,">>$userlog" or die("$userlog PATH ERROR!\n");
		print OUT1 "$detail \n";
		close OUT1;
	}else{
		print "$detail \n";
	}
}
sub mkdirOrDie{
	my ($dir) = @_ ;
	if(!-d $dir){
		make_path($dir);
		$dir=abs_path($dir);
	}
}
