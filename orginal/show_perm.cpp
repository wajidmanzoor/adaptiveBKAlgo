#include "inc/common.h"
#include "inc/graph.h"

int main(int argc, const char *argv[]) {
  string path(argv[1]);
  Graph g(path);
  ui n = g.n;
  vector<ui> deg(g.degree.begin(), g.degree.end());
  ui maxDeg = *max_element(deg.begin(), deg.end());

  vector<ui> bins(maxDeg+1,0);
  for (ui d : deg) bins[d]++;
  vector<ui> binStart(maxDeg+1,0);
  partial_sum(bins.begin(), bins.end()-1, binStart.begin()+1);

  vector<ui> pos(n), sorted(n);
  for (ui v=0;v<n;v++) { pos[v]=binStart[deg[v]]++; sorted[pos[v]]=v; }
  for (ui d=0;d<=maxDeg;d++) binStart[d]-=bins[d];

  vector<ui> peelSeq(n);
  vector<ui> core(n);
  for (ui i=0;i<n;i++) {
    ui v=sorted[i];
    core[v]=deg[v];
    peelSeq[n-1-i]=v;
    for (ui j=g.offset[v];j<g.offset[v+1];j++) {
      ui u=g.neighbors[j];
      if (deg[u]>deg[v]) {
        ui du=deg[u],pu=pos[u],pw=binStart[du],w=sorted[pw];
        if (u!=w) { pos[u]=pw;sorted[pu]=w;pos[w]=pu;sorted[pw]=u; }
        binStart[du]++; deg[u]--;
      }
    }
  }

  // ASC: peelSeq[n-1-i] gets new label i (first peeled = label 0)
  vector<ui> perm(n);
  for (ui i=0;i<n;i++) perm[peelSeq[n-1-i]]=i;

  printf("orig_id -> new_label  (core_value, orig_degree)\n");
  for (ui v=0;v<n;v++)
    printf("  orig %2d  ->  new %2d   (core=%d, deg=%d)\n",
           v, perm[v], core[v], (int)g.degree[v]);

  printf("\nnew_label -> orig  with expandTo in new labels\n");
  // collect per new label
  vector<ui> origOf(n);
  for (ui v=0;v<n;v++) origOf[perm[v]]=v;

  for (ui newv=0;newv<n;newv++) {
    ui v=origOf[newv];
    vector<ui> higher;
    for (ui j=g.offset[v];j<g.offset[v+1];j++) {
      ui u=g.neighbors[j];
      if (perm[u]>newv) higher.push_back(perm[u]);
    }
    sort(higher.begin(),higher.end());
    printf("  new %2d (orig %2d, core=%d):  expandTo={ ", newv, v, core[v]);
    for (ui x:higher) printf("%d ",x);
    printf("}\n");
  }
  return 0;
}
