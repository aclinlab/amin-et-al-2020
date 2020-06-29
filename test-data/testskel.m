s = neurSkel('testskel.csv');
s = s.getSynSets({'testsyn.csv', 'testsyn-2.csv'}, [1 2]);
s = s.getSynSetGraph(1);
s = s.getSynSetGraph(2);
