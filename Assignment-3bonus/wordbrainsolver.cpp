// Copyright year Negmetulla_Yerlan 180107219@stu.sdu.edu.kz
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <stdlib.h>
using namespace std;

const short int LETTERS_SIZE = 26;


struct Trie {
    char path;
    struct Trie *path_0;
    struct Trie *letter[LETTERS_SIZE];
    string word = "";
    bool isEndOfWord = false;
};


struct wordlist {
    int n = 0;
    vector<vector<char>> wordlist_0;
    void wordlist_1(int m) {
        n = m;
        wordlist_0.clear();
        wordlist_0.resize(n, vector<char> (m, 0));
    }
};


bool checkT(
    vector<string> tree, 
    wordlist *point,
    vector<int> len,
    Trie *word_trie,
    map<int, char> *Map,
    int word_file0, int word_file1,
    vector<string> *words,
    set<string> *map_set
);


bool checkA(
    vector<string> tree, 
    wordlist *point,
    Trie *word_trie,
    vector<int> len,
    map<int, char> Map,
    vector<string> *words,
    set<string> *map_set
){
    if (len.size() == 0) {
        string str_tmp = "";
        for (const auto &pointer : *words)
        str_tmp += pointer + " ";
        map_set->insert(str_tmp);
        return 0;
    }
    wordlist word_file_trie = *point;
    for (
        map<int, char>::reverse_iterator Map_new = Map.rbegin();
        Map_new != Map.rend(); 
        ++Map_new
    ) {
        word_file_trie.wordlist_0[Map_new->first % word_file_trie.n].erase(
            word_file_trie.wordlist_0[Map_new->first % word_file_trie.n].begin() 
            + Map_new->first / word_file_trie.n
        );
    }
    for (int i = 0; i < word_file_trie.wordlist_0.size(); i++) {
        for (int j = 0; j < word_file_trie.wordlist_0.at(i).size(); j++) {
            map<int, char> map_trie;
            vector<string> words_trie = *words;
            string str_0, str_1;
            if (words_trie.size() != 0) {
                for (int x = 0; x < words_trie.size(); x++) {
                    str_0 = words_trie[x];
                    str_1 = tree[x];
                    int count = 0;
                    for (int y = 0; y < str_0.size(); y++) {
                        if (str_1[y] != '*') {
                            if (str_0[y] == str_1[y])
                                count++;
                        } else {
                            count++;
                        }
                    }
                    if (count != str_0.size())
                        return 1;
                }
            }
            if (word_trie->letter[len.front()]) {
                if (checkT(
                    tree, 
                    &word_file_trie, 
                    len,
                    word_trie->letter[len.front()]->letter[word_file_trie.
                    wordlist_0.at(i).at(j) - 'a'],
                    &map_trie, 
                    i, j, 
                    &words_trie, 
                    map_set
                )) {
                    *words = words_trie;
                    return 1;
                }
            }
        }
    }
    return 0;
}
Trie *word_trie_0;


bool checkT(
    vector<string> tree, 
    wordlist *point,
    vector<int> len,
    Trie *word_trie,
    map<int, char> *Map,
    int word_file0, int word_file1,
    vector<string> *words,
    set<string> *map_set
){
    if (!word_trie) return 0;
    map<int, char> map_trie;
    map_trie.insert(Map->begin(), Map->end());
    map_trie[word_file0 + word_file1 * point->n] = word_trie->path;
    for (const auto &pointer : (*Map))
        if (word_trie->isEndOfWord) {
            vector<int> trie_len = len;
            trie_len.erase(trie_len.begin());
            vector<string> words_trie = *words;
            words_trie.push_back(word_trie->word);
            for (const auto &pointer : (words_trie))
                if (checkA(tree, point, word_trie_0, trie_len, map_trie, &words_trie, map_set)) {
                    *words = words_trie;
                    *Map = map_trie;
                    return 1;
                } else return 0;
        }
    int begin_0 = (word_file0 - 1 > 0) ? word_file0 - 1 : 0;
    int end_0 = (word_file1 - 1 > 0) ? word_file1 - 1 : 0;
    int begin_1 = (word_file0 + 1 < point->n) ? word_file0 + 1 : point->n - 1;
    int end_1 = (word_file1 + 1 < point->n) ? word_file1 + 1 : point->n - 1;
    for (int i = begin_0; i <= begin_1; i++)
        for (int j = end_0; j <= end_1; j++) {
            if (map_trie.count(i + j * point->n)) continue;
            char wt;
            try {
                wt = point->wordlist_0.at(i).at(j);
            } catch(out_of_range) {
                continue;
            }
            if (word_trie->letter[wt - 'a'])
                if (checkT(
                    tree, 
                    point, 
                    len, 
                    word_trie->letter[point->wordlist_0.at(i).at(j) - 'a'], 
                    &map_trie, 
                    i, j, 
                    words, 
                    map_set
                )) {
                    Map->insert(map_trie.begin(), map_trie.end());
                    return 1;
                }
        }
    return 0;
}


int read_file(string ss, Trie *word_trie) {
    ifstream infile(ss);
    string word_file;
    int count;
    while (infile >> word_file) {
        if (word_file.length() > LETTERS_SIZE) continue;
        count++;
        if (word_trie->letter[word_file.length()] == nullptr) {
            word_trie->letter[word_file.length()] = new Trie;
            word_trie->letter[word_file.length()]->path_0 = word_trie;
            word_trie->letter[word_file.length()]->path = '0' + word_file.length();
        }
        Trie *bypasses = word_trie->letter[word_file.length()];
        for (auto wt : word_file) {
            if (bypasses->letter[wt-'a'] == nullptr) {
                bypasses->letter[wt-'a'] = new Trie;
                bypasses->letter[wt-'a']->path_0 = bypasses;
                bypasses->letter[wt-'a']->path = wt;
            }
            bypasses = bypasses->letter[wt - 'a'];
        }
        bypasses->isEndOfWord = true;
        bypasses->word = word_file;
    }
    return count;
}


int main(int argc, char **argv) {
    Trie word_trie = Trie();
    Trie word_trie_len = Trie();
    read_file(argv[1], &word_trie);
    read_file(argv[2], &word_trie_len);
    while (true) {
        int flag = 0;
        wordlist point = wordlist();
        string str;
        vector<vector<char>> result;
        string word_file;
        string line = "";
        while (cin) {
            getline(cin, line);
            if (line == "") exit(0);
            if (line.find("*") == string::npos) {
                vector<char> data(line.begin(), line.end());
                result.push_back(data);
            } else {
                word_file = line;
                break;
            }
        }
        int tree_len = result[0].size();

        vector<string> tree;
        string::size_type tmp, mark;
        mark = word_file.find(" ");
        tmp = 0;
        while (string::npos != mark) {
            tree.push_back(word_file.substr(tmp, mark - tmp));
            tmp = mark + string(" ").size();
            mark = word_file.find(" ", tmp);
        }
        if (tmp != word_file.length())
            tree.push_back(word_file.substr(tmp));

        word_trie_0 = &word_trie;
        map<int, char> Map;
        vector<int> len;
        
        point.wordlist_1(tree_len);
        for (int i = 0; i < point.n; i++) {
            for (int j = 0; j < result.size(); j++)
                point.wordlist_0[i][j] = result[result.size() - 1 - j][i];
        }
        for (auto &k : tree)
            len.push_back(k.length());
        if (len.size() == 0 || point.n == 0)
            return 0;

        vector<string> words;
        set<string> map_set;
        checkA(tree, &point, &word_trie, len, Map, &words, &map_set);
        if (map_set.size() == 0) {
            word_trie_0 = &word_trie_len;
            checkA(
                tree,
                &point, 
                &word_trie_len, 
                len, 
                Map, 
                &words, 
                &map_set
            );
        }
        for (const auto &pointer : map_set)
            cout << pointer << endl;
        cout << "." << endl;
    }
}
