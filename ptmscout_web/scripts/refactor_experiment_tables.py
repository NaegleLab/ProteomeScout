from paste.deploy.loadwsgi import appconfig
import os
import traceback
from ptmscout.database import experiment, DBSession
from DB_init import DatabaseInitialization
import re


def get_month(m):
    months = ['january','february','march','april','may','june','july','august','september','october', 'november','december']
    
    if(m == ''):
        return m
    
    try:
        i = int(m)
        return months[i-1]
    except:
        for month in months:
            if month.find(m.lower())==0:
                return month
    return None

if __name__ == '__main__':
    try:
        settings = appconfig(os.path.join('config:', 'data', 'ptmscout', 'ptmscout_web', 'development.ini'))
            
        DatabaseInitialization.setUpClass(settings)
        dbinit = DatabaseInitialization()
        dbinit.setUp()
        
        
        experiments = DBSession.query(experiment.Experiment).all()
        
        for exp in experiments:
            pg = exp.pages
            pub_date = exp.pub_date
            
            try:
                regex = re.compile("(?P<section>[A-Za-z]+)?(?P<start>[0-9]+)(?:\-(?P<end>[A-Za-z]*[0-9]+))?")
                matcher = regex.match(pg)
                
                gdict = matcher.groupdict()
                
                section = None
                start = None
                end = None
                
                if('section' in gdict):
                    section = gdict['section']
                if('start' in gdict):
                    start = gdict['start']
                if('end' in gdict):
                    end = gdict['end']
                
                if start == None:
                    raise Exception("No start page")
                
                if end != None and int(end) < int(start):
                    header = start[0:-len(end)]
                    end = header + end
                
                if(section == None): section = ""
                
                start = section+start
                if end != None:
                    end = section + end
                
                exp.page_start = start
                exp.page_end = end
                
                print "Pages:", pg, start, end
                
            except Exception, e:
                print "Couldn't parse page limits: " + pg, str(e)

            try:
                regex_text = re.compile("([0-9]{4})\-([A-Za-z]+)")
                regex_num = re.compile("([0-9]{4})(?:\-([0-9]{2}))?")
                
                matcher_text = regex_text.match(pub_date)
                matcher_num = regex_num.match(pub_date)
                
                if(matcher_text != None):
                    year = int(matcher_text.group(1))
                    month = matcher_text.group(2)
                elif(matcher_num != None):
                    year = int(matcher_num.group(1))
                    month = matcher_num.group(2)
                    if(month == None):
                        month = ''
                else:
                    raise Exception()
                
                m = get_month(month)
                if m == None:
                    raise Exception("No such month: " + str(month))
                
                exp.publication_year = year
                exp.publication_month = m
                
                print " Date:", pub_date, year, m
                
            except Exception, e:
                print "Couldn't parse date: " + pub_date, str(e)
        DBSession.flush()
    except Exception, e:
        traceback.print_exc()
        
        print "Rolling back database changes..."
        dbinit.rollback()
    else:
        print "Finalizing DB changes"
        dbinit.tearDown()