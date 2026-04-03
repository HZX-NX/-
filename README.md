#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>
#include <random>
#include <map>

using namespace std;

// ---------- 基因数据库 ----------
struct Gene {
    string pattern;
    string name;
    string species;
    string function;
};

vector<Gene> geneDatabase = {
    {"ATGGTGCTGTCCGCTGCCGAC", "HBB (β-globin)", "Homo sapiens", "编码血红蛋白β亚基，参与氧气运输。"},
    {"ATGGCGCTGCTCACAGAGTTC", "TP53", "Homo sapiens", "肿瘤抑制蛋白p53，调控细胞周期和凋亡。"},
    {"ATGGCCAAGGTGAAGGTCGGAGTC", "GAPDH", "Homo sapiens", "甘油醛-3-磷酸脱氢酶，糖酵解关键酶。"},
    {"ATGGAAATAGCTTGTGAAAAGAC", "BRCA2", "Homo sapiens", "乳腺癌易感蛋白2，参与同源重组修复。"},
    {"ATGCAGAGGTCGCCTCTGGAAAAGGCCAGC", "CFTR", "Homo sapiens", "囊性纤维化跨膜传导调节因子，氯离子通道。"},
    {"ATGGAAGAAAAGAGTCGCCGCA", "EGFR", "Homo sapiens", "表皮生长因子受体，受体酪氨酸激酶。"},
    {"ATGACTGAATATAAACTTGTGGTAGTTGGAGCT", "KRAS", "Homo sapiens", "GTP酶，参与细胞信号转导。"},
    {"ATGGACTTTTCGCAGGTG", "MYC", "Homo sapiens", "转录因子，调控细胞增殖和凋亡。"}
};

// ---------- 氨基酸密码子表 ----------
map<string, char> codonMap = {
    {"ATA",'I'},{"ATC",'I'},{"ATT",'I'},{"ATG",'M'},
    {"ACA",'T'},{"ACC",'T'},{"ACG",'T'},{"ACT",'T'},
    {"AAC",'N'},{"AAT",'N'},{"AAA",'K'},{"AAG",'K'},
    {"AGC",'S'},{"AGT",'S'},{"AGA",'R'},{"AGG",'R'},
    {"CTA",'L'},{"CTC",'L'},{"CTG",'L'},{"CTT",'L'},
    {"CCA",'P'},{"CCC",'P'},{"CCG",'P'},{"CCT",'P'},
    {"CAC",'H'},{"CAT",'H'},{"CAA",'Q'},{"CAG",'Q'},
    {"CGA",'R'},{"CGC",'R'},{"CGG",'R'},{"CGT",'R'},
    {"GTA",'V'},{"GTC",'V'},{"GTG",'V'},{"GTT",'V'},
    {"GCA",'A'},{"GCC",'A'},{"GCG",'A'},{"GCT",'A'},
    {"GAC",'D'},{"GAT",'D'},{"GAA",'E'},{"GAG",'E'},
    {"GGA",'G'},{"GGC",'G'},{"GGG",'G'},{"GGT",'G'},
    {"TCA",'S'},{"TCC",'S'},{"TCG",'S'},{"TCT",'S'},
    {"TTC",'F'},{"TTT",'F'},{"TTA",'L'},{"TTG",'L'},
    {"TAC",'Y'},{"TAT",'Y'},{"TAA",'*'},{"TAG",'*'},
    {"TGC",'C'},{"TGT",'C'},{"TGA",'*'},{"TGG",'W'}
};

// ---------- 辅助函数 ----------
string toUpper(string s) {
    transform(s.begin(), s.end(), s.begin(), ::toupper);
    return s;
}

string reverseComplement(const string& seq) {
    map<char, char> comp = {{'A','T'},{'T','A'},{'C','G'},{'G','C'}};
    string rc = seq;
    reverse(rc.begin(), rc.end());
    for (char& c : rc) c = comp[c];
    return rc;
}

double gcContent(const string& seq) {
    int gc = count(seq.begin(), seq.end(), 'G') + count(seq.begin(), seq.end(), 'C');
    return (double)gc / seq.length() * 100.0;
}

// 查找开放阅读框 (最长 ORF)
pair<int, string> findLongestORF(const string& seq) {
    int bestStart = -1, bestEnd = -1, bestLen = 0;
    for (size_t i = 0; i + 2 < seq.length(); ++i) {
        if (seq.substr(i, 3) == "ATG") {
            for (size_t j = i+3; j + 2 < seq.length(); j += 3) {
                string codon = seq.substr(j, 3);
                if (codon == "TAA" || codon == "TAG" || codon == "TGA") {
                    int len = j - i + 3;
                    if (len > bestLen) {
                        bestLen = len;
                        bestStart = i;
                        bestEnd = j+2;
                    }
                    break;
                }
            }
        }
    }
    if (bestStart != -1)
        return {bestStart, seq.substr(bestStart, bestLen)};
    else
        return {-1, ""};
}

// 翻译 ORF 为蛋白质序列
string translate(const string& orf) {
    string protein;
    for (size_t i = 0; i + 2 < orf.length(); i += 3) {
        string codon = orf.substr(i, 3);
        if (codonMap.count(codon))
            protein += codonMap[codon];
        else
            protein += '?';
    }
    return protein;
}

// 匹配基因数据库
const Gene* findMatchingGene(const string& seq) {
    for (const auto& gene : geneDatabase) {
        if (seq.find(gene.pattern) != string::npos)
            return &gene;
    }
    return nullptr;
}

// 模拟 Sanger 测序曲线图（文本版）
void printSangerTrace(const string& seq) {
    cout << "\n--- Sanger 测序四色荧光图谱（文本模拟）---\n";
    int len = seq.length();
    for (int i = 0; i < len; ++i) {
        cout << seq[i];
        if ((i+1) % 60 == 0) cout << "\n";
        else if ((i+1) % 10 == 0) cout << " ";
    }
    cout << "\n每个碱基对应的峰高（随机模拟，仅供示意）：\n";
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(30, 100);
    for (int i = 0; i < len; ++i) {
        int height = dis(gen);
        if (i > 0 && seq[i] == seq[i-1]) height *= 0.8;
        cout << seq[i] << ":" << height << "  ";
        if ((i+1) % 20 == 0) cout << "\n";
    }
    cout << "\n";
}

// 生成 GenBank 风格注释
void printGenBankAnnotation(const string& seq) {
    const Gene* gene = findMatchingGene(seq);
    if (gene) {
        cout << "LOCUS       " << gene->name << "                " << seq.length() << " bp    DNA    linear   "
             << (gene->species == "Homo sapiens" ? "HUM" : "UNK") << "\n";
        cout << "DEFINITION  " << gene->species << " " << gene->name << " gene, partial cds.\n";
        cout << "ACCESSION   SIM_" << gene->name << "_1\n";
        cout << "VERSION     SIM.1\n";
        cout << "SOURCE      " << gene->species << "\n";
        cout << "  ORGANISM  " << gene->species << "\n";
        cout << "FEATURES             Location/Qualifiers\n";
        cout << "     CDS             1.." << seq.length() << "\n";
        cout << "                     /gene=\"" << gene->name << "\"\n";
        cout << "                     /product=\"" << gene->function << "\"\n";
        cout << "ORIGIN\n";
        for (size_t i = 0; i < seq.length(); i += 60) {
            cout << seq.substr(i, 60) << "\n";
        }
        cout << "//\n\n";
        cout << "📌 功能描述：" << gene->function << "\n";
        cout << "🔬 物种：" << gene->species << "\n";
        cout << "📖 基因名称：" << gene->name << "\n";
    } else {
        cout << "LOCUS       UNKNOWN                 " << seq.length() << " bp    DNA    linear   UNK\n";
        cout << "DEFINITION  Uncharacterized sequence.\n";
        cout << "FEATURES             Location/Qualifiers\n";
        cout << "     misc_feature    1.." << seq.length() << "\n";
        cout << "                     /note=\"Unknown gene\"\n";
        cout << "ORIGIN\n";
        for (size_t i = 0; i < seq.length(); i += 60) {
            cout << seq.substr(i, 60) << "\n";
        }
        cout << "//\n\n";
        cout << "📌 该序列未匹配到已知基因。\n";
    }
}

// ---------- 主程序 ----------
int main() {
    cout << "========================================\n";
    cout << "   基因测序模拟器 (C++ 命令行版)\n";
    cout << "========================================\n";
    cout << "输入 DNA 序列 (5'→3')，仅含 A/T/C/G：\n";
    string seq;
    getline(cin, seq);
    seq = toUpper(seq);
    // 移除空格
    seq.erase(remove(seq.begin(), seq.end(), ' '), seq.end());
    if (seq.find_first_not_of("ATCG") != string::npos) {
        cout << "错误：序列只能包含 A, T, C, G 字符。\n";
        return 1;
    }

    // 1. GC 含量
    double gc = gcContent(seq);
    cout << "\n🧬 GC 含量: " << gc << "%\n";
    cout << "📏 序列长度: " << seq.length() << " bp\n";
    cout << "🔄 反向互补: " << reverseComplement(seq) << "\n";

    // 2. Sanger 图谱文本模拟
    printSangerTrace(seq);

    // 3. ORF 查找与翻译
    auto [orfStart, orfSeq] = findLongestORF(seq);
    if (orfStart != -1) {
        string protein = translate(orfSeq);
        cout << "\n🧬 **开放阅读框 (ORF)**\n";
        cout << "起始于第 " << orfStart+1 << " 位，终止于第 " << orfStart+orfSeq.length() << " 位\n";
        cout << "ORF 序列: " << orfSeq << "\n";
        cout << "翻译蛋白 (5'→3'): " << protein << "\n";
    } else {
        cout << "\n🧬 未发现有效 ORF（起始密码子后无终止密码子）\n";
    }

    // 4. GenBank 风格注释
    cout << "\n📖 基因解读 (NCBI/GenBank 风格注释)\n";
    printGenBankAnnotation(seq);

    return 0;
}
