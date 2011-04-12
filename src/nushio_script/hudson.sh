rm -fr ~/comprehensive_test
ssh t2a006180 'mkdir -p ~/comprehensive_test'
cd ~/comprehensive_test/
svn co svn://130.54.55.241/home/svn/repos/nushiolab2/mhd/2010-12 .
rake link
ssh t2a006180 'cd ~/comprehensive_test ;rake test'
ssh t2a006180 'cd ~/comprehensive_test ;rake'



