# 螺旋状展開（Spiral Unfolding）に関する既往研究レビュー

---

## 1. はじめに

多面体や多胞体の**展開図（net / unfolding）**を構成する問題は，
15世紀のデューラーに遡る古典的な問題であると同時に，
計算幾何学・組合せ論における現代的な研究課題でもある．
本レビューでは，Apple-Peel 展開に関連する既往研究を以下の5つの観点から整理する：
(1) デューラー問題と基礎的展開手法，
(2) ハミルトン的・Zipper 展開，
(3) 螺旋展開・Apple-Peel 展開，
(4) Band 展開，
(5) 4次元多胞体の展開．

---

## 2. デューラー問題と基礎的展開手法

### 2.1 問題の起源と定式化

多面体の展開図の起源はルネサンス期の芸術家 **Albrecht Dürer** の著作
*Underweysung der Messung* (1525) に求められる．
この著作には立方体や正十二面体などの展開図の手描き図が含まれており，
「辺に沿って切ることで重なりのない展開図が必ずできるか」という問いの萌芽が見られる．
この問いが数学的命題として定式化されたのは約450年後，
**Shephard (1975)** によってである．

> **Shephard (1975)**
> *Convex Polytopes with Convex Nets*
> Mathematical Proceedings of the Cambridge Philosophical Society, Vol. 78, pp. 389–403
> DOI: 10.1017/S0305004100051860

Shephard は「すべての凸多面体は辺に沿った切断によって重なりのない展開図を持つか」を
命題として提示した（現在も未解決）．
また，Hamiltonian な辺パスに沿った切断による展開（Hamiltonian unfolding）の概念も導入した．

### 2.2 保証された展開手法

凸多面体に対して重なりなしを**保証**する展開手法として，以下の2つが知られている．

> **Aronov & O'Rourke (1992)**
> *Nonoverlap of the Star Unfolding*
> Discrete & Computational Geometry, Vol. 8, No. 3, pp. 219–250
> DOI: 10.1007/BF02293047

**Star unfolding** は，ある源点から各頂点への測地最短路に沿って切断する手法であり，
任意の凸多面体に対して重なりのない展開図が得られることが証明された．
Source unfolding（Voronoi 図に基づく切断）とともに，凸多面体に対する
2つの保証付き展開手法として位置づけられる．

### 2.3 アフィン変換による進展

> **Ghomi (2014)**
> *Affine Unfoldings of Convex Polyhedra*
> Geometry & Topology, Vol. 18, No. 5, pp. 3055–3090
> arXiv:1305.3231

Ghomi は，一般の位置にある凸多面体はアフィン変換（線形引き伸ばし）を施した後に
辺展開可能であることを証明した．
これはデューラー問題に組合せ論的障害は存在しないことを示す最も重要な近年の進展である．

### 2.4 包括的参考書

> **Demaine & O'Rourke (2007)**
> *Geometric Folding Algorithms: Linkages, Origami, Polyhedra*
> Cambridge University Press

展開・折り畳みに関する計算幾何学の決定版教科書．
展開図・デューラー問題・Cauchy の剛性定理・Alexandrov の一意性定理など
関連するすべての話題を網羅し，多数の未解決問題を列挙する．

---

## 3. ハミルトン的展開・Zipper 展開

「螺旋状に展開する」というアイデアに最も近いのが，
**Hamiltonian unfolding（Zipper unfolding）** の系統である．

> **Lubiw, Demaine et al. (2010)**
> *Zipper Unfoldings of Polyhedral Complexes*
> Proceedings of CCCG 2010, pp. 219–222

**Zipper unfolding** とは，すべての面を1本の連続するハミルトン辺パスで切断して
得られる展開図である（Apple-Peel の辺切断版に相当）．
すべての Platonic 立体および Archimedean 立体が Hamiltonian 展開（Zipper 展開）を
持つことが示されており，Apple-Peel 展開の直接の先行概念といえる．

> **Demaine, Demaine & Uehara (2013)**
> *Zipper Unfolding of Domes and Prismoids*
> Proceedings of CCCG 2013

Dome 型多面体に対して，グラフがハミルトン的であっても
ハミルトン展開がすべて重なりを生じる例を示す一方，
nested prismoid はハミルトン展開可能であることを証明した．

> **O'Rourke (2010)**
> *Flat Zipper-Unfolding Pairs for Platonic Solids*
> arXiv:1010.2450; CCCG 2010

正四面体・立方体・正八面体・正二十面体の4種は，
ハミルトン辺パスで切断した展開図が平行四辺形の二重被覆に折り畳める
「Zipper ペア」を持つことを示す．正十二面体では構成不可能．

---

## 4. 螺旋展開・Apple-Peel 展開

### 4.1 Spiral Unfolding（O'Rourke, 2015）

Apple-Peel 展開に直接対応する先行研究として，O'Rourke による以下の論文がある．

> **O'Rourke (2015)**
> *Spiral Unfoldings of Convex Polyhedra*
> arXiv:1509.00321

本論文は，凸多面体のすべての面を螺旋状の帯で覆う切断パス
（**spiral unfolding**）を定義し，その重なりの有無を体系的に研究した．
Platonic・Archimedean 立体13種すべてが重なりのない螺旋展開を持つことを証明し，
回転体（polyhedra of revolution）を主要な解析クラスとした．
一般の凸多面体では螺旋展開が重なりを持つ場合があることも示している．
本プロジェクトの Apple-Peel 展開は，面の選択を連続的な帯ではなく
**離散的なグリーディ選択**に置き換えた発展版と位置づけられる．

### 4.2 Apple-Peel Unfolding（Yoshino & Chaidee）

> **Yoshino & Chaidee (~2024)**
> *Apple Peel Unfolding of Archimedean and Catalan Solids*
> Submitted to Journal of Mathematics and Arts

Archimedean 立体13種とその双対（Catalan 立体）に対して Apple-Peel 展開を構成し，
重なりなし条件を検討する論文（投稿中）．
本プロジェクトの3次元アルゴリズムの直接的な出典である．

### 4.3 日本語圏における先行研究

> **Unno et al. (2012)**
> *リンゴの皮むきと正多面体の展開図について*
> 形の科学会誌, Vol. 27, No. 1, pp. 27–28

りんごを剥くときの螺旋パターンと正多面体の展開図の幾何学的関係を論じた
日本語初期論文．Apple-Peel の発想との共鳴が見られる．

---

## 5. Band 展開

螺旋展開と関連する別の系統として，**band unfolding** がある．
多面体の表面のうち2つの平行な平面に挟まれた帯状の部分（polyhedral band）の展開を研究する．

> **Aloupis, Demaine et al. (2008)**
> *Edge-Unfolding Nested Polyhedral Bands*
> Computational Geometry: Theory and Applications, Vol. 39, No. 1, pp. 30–42

Nested band（一方の境界多角形が他方を正射影に含む構造）は
常に辺展開可能であることを証明した基礎的結果．
螺旋展開における「帯を一周ずつ展開する」動作の理論的根拠に対応する．

> **Radons (2023)**
> *Edge-Unfolding Nested Prismatoids*
> Computational Geometry: Theory and Applications

Nested prismatoid（2つの平行凸多角形の凸包）は必ず辺展開可能であることを証明し，
O'Rourke (2007) の反例で残された未解決問題を解決した．

---

## 6. 4次元多胞体の展開

### 6.1 組合せ論的計数

> **Buekenhout & Parker (1998)**
> *The Number of Nets of the Regular Convex Polytopes in Dimension ≤ 4*
> Discrete Mathematics, Vol. 186, pp. 69–94

次元 ≤ 4 のすべての正凸多胞体について，
面隣接グラフの全域木（対称性の同値類）を分類・計数した基礎的論文．
最も著名な結果：**超立方体（tesseract）の展開図は261種**（Turney 1984 の計算を厳密化）．

### 6.2 4次元正多胞体の展開可能性

> **Devadoss & Harvey (2022)**
> *Unfoldings and Nets of Regular Polytopes*
> Computational Geometry: Theory and Applications
> arXiv:2111.01359

$n$-次元立方体・$n$-単体・4-直交胞（4-orthoplex）は
「all-net」（すべての稜展開が有効な展開図を与える）であることを証明した．
一方，5次元以上の直交胞や600胞体では反例（重なりを持つ展開図）が存在することを示す．
本プロジェクトが対象とする多胞体の展開可能性理論における重要な参照点である．

> **Devadoss, Harvey & Zhang (2022)**
> *Visualizing and Unfolding Nets of 4-Polytopes*
> SoCG 2022 (Media Exposition), LIPIcs Vol. 224, Article 67

4次元多胞体の展開図を対話的に可視化するウェブツールを発表した Media Exposition 論文．
4-単体・4-立方体・4-直交胞を対象とし，双対グラフ上の全域木を描くことで
展開図を段階的に構築できる．

### 6.3 4次元における Hamiltonian 展開

> **Akitaya & Samanta et al. (2024)**
> *Path-Unfolding the Tesseract*
> Fall Workshop on Computational Geometry (FWCG 2024)

Tesseract（超立方体）の**path unfolding**（双対グラフがハミルトンパスになる展開）を
すべて列挙した論文（合計 35,520 通り）．
3次元の Zipper/Hamiltonian 展開を4次元に直接拡張したものであり，
本プロジェクトの Apple-Peel 展開の4次元版との概念的一致が大きい．

### 6.4 4次元多胞体のハミルトン的構造

> **Séquin (2005)**
> *Symmetrical Hamiltonian Manifolds on Regular 3D and 4D Polytopes*
> Bridges Conference 2005, pp. 463–472

正多胞体の辺グラフ上のハミルトン的構造（サイクルから多様体への拡張）を研究した論文．
4-単体ではメビウス帯状の三角形被覆が得られることなどを示す．
螺旋展開と4次元のハミルトン構造の接点を考察するうえで参考になる．

---

## 7. 考察：Apple-Peel 展開の学術的位置づけ

以上の既往研究を踏まえると，Apple-Peel 展開は以下のように位置づけられる．

| 観点 | 既往研究 | Apple-Peel 展開 |
|:---|:---|:---|
| 切断の単純性 | 辺切断のみ（Zipper），面切断も可（一般） | 面の選択順序（切断線は明示しない） |
| 螺旋性の保証 | O'Rourke (2015)：Platonic/Archimedean で保証 | グリーディ選択で螺旋性を実現（保証は一部） |
| 対称性への依存 | Band 展開：nested 構造が必要 | 対称性が高いほど成功率が高い（Archimedean で実証） |
| 4次元への拡張 | Devadoss (2022)：理論，Akitaya (2024)：列挙 | 本プロジェクト：グリーディ + Darboux フレーム |
| アルゴリズム | Tabu 探索（Zawallich 2024）など | グリーディ + 1段バックトラック |

特に重要な相違点は，Apple-Peel 展開が「どの面を次に選ぶか」という
**選択規則（greedy rule）** を明示的に定義している点であり，
Darboux フレームによる左方向の定式化（本プロジェクトで導入した行列式条件）は
従来の螺旋展開研究には見られない独自の貢献である．

また，**4次元への拡張**については，Devadoss & Harvey (2022) が全域木による
一般的な展開を論じているのに対し，本プロジェクトは螺旋的な**選択順序**に焦点を当てており，
Akitaya et al. (2024) の path-unfolding の精神とより近い．
Apple-Peel 展開の4次元版は，単なる展開図の列挙ではなく，
「螺旋」の幾何学的意味を4次元で定式化するという新規な試みである．

---

## 参考文献一覧

| 文献 | 著者 | 年 | 掲載誌／媒体 |
|:---|:---|:---|:---|
| *Underweysung der Messung* | Dürer | 1525 | 著書 |
| Convex Polytopes with Convex Nets | Shephard | 1975 | Math. Proc. Cambridge |
| Nonoverlap of the Star Unfolding | Aronov, O'Rourke | 1992 | Discrete & Comput. Geom. |
| The Number of Nets of Regular Convex Polytopes in Dim ≤ 4 | Buekenhout, Parker | 1998 | Discrete Mathematics |
| Unfolding Some Classes of Orthogonal Polyhedra | Biedl, Demaine et al. | 1998 | CCCG 1998 |
| Ununfoldable Polyhedra with Convex Faces | Bern, Demaine et al. | 2003 | Comput. Geom.: Theory Appl. |
| Geometric Folding Algorithms (書籍) | Demaine, O'Rourke | 2007 | Cambridge Univ. Press |
| Edge-Unfolding Nested Polyhedral Bands | Aloupis, Demaine et al. | 2008 | Comput. Geom.: Theory Appl. |
| Symmetrical Hamiltonian Manifolds on 3D and 4D Polytopes | Séquin | 2005 | Bridges 2005 |
| Zipper Unfoldings of Polyhedral Complexes | Lubiw, Demaine et al. | 2010 | CCCG 2010 |
| Flat Zipper-Unfolding Pairs for Platonic Solids | O'Rourke | 2010 | arXiv:1010.2450 |
| リンゴの皮むきと正多面体の展開図 | Unno et al. | 2012 | 形の科学会誌 |
| Zipper Unfolding of Domes and Prismoids | Demaine, Demaine, Uehara | 2013 | CCCG 2013 |
| Affine Unfoldings of Convex Polyhedra | Ghomi | 2014 | Geom. & Topology |
| Spiral Unfoldings of Convex Polyhedra | O'Rourke | 2015 | arXiv:1509.00321 |
| A Survey of Folding and Unfolding | Demaine, O'Rourke | 2005/2007 | MSRI Publ. Vol. 52 |
| Dürer's Unfolding Problem for Convex Polyhedra | Ghomi | 2018 | Notices AMS |
| Unfolding Polyhedra (survey) | O'Rourke | 2019 | CCCG 2019 / arXiv |
| Apple Peel Unfolding of Archimedean and Catalan Solids | Yoshino, Chaidee | ~2024 | J. Math. Arts (投稿中) |
| Unfoldings and Nets of Regular Polytopes | Devadoss, Harvey | 2022 | Comput. Geom.: Theory Appl. |
| Visualizing and Unfolding Nets of 4-Polytopes | Devadoss, Harvey, Zhang | 2022 | SoCG 2022 |
| Edge-Unfolding Nested Prismatoids | Radons | 2023 | Comput. Geom.: Theory Appl. |
| Unfolding Polyhedra via Tabu Search | Zawallich | 2024 | The Visual Computer |
| Path-Unfolding the Tesseract | Akitaya, Samanta et al. | 2024 | FWCG 2024 |

---

*本レビューにおいて，未発表・投稿中の文献（Yoshino & Chaidee ~2024）については
情報が限定的であり，詳細は著者に確認が必要である．
また，Kaino (2019) の4次元展開図に関する日本語発表（対称性学会 金沢大会）は，
国際的なデータベースへの収録が確認できなかったため本レビューへの掲載を見送った．*
