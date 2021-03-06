(ns rename-genes)

(import '(java.util.regex Matcher))


(defn is-comment-line? [xs]
 (=
   (first (first xs))
   \#))

(def new-id (atom 0))

(defn gen-new-id [anno species-abbrev]
  (str species-abbrev "." anno "." (swap! new-id inc)))

(defn rename-gene [line prev-id new-id]
 (clojure.string/replace line (Matcher/quoteReplacement prev-id) new-id))

(defn annotate-first-line [line]
  (str line ";annotation_note=\"Generated by a second round of gene prediction using AUGUSTUS-PPX\""))

(with-open [rdr (clojure.java.io/reader "predicted_genes.unprocessed.gff3")]
  (let [split-lines (partition-by (fn [x] (re-find #"^#" x)) (line-seq rdr))
        [family-name species-abbrev] *command-line-args*]
    (doseq [lines (remove is-comment-line? split-lines)]
      (let [gene-id (second (re-find #"ID=(.+)" (first lines)))
            new-id (gen-new-id family-name species-abbrev)
            first-line (annotate-first-line (rename-gene (first lines) gene-id new-id))
            rest-lines (map (fn [x] (rename-gene x gene-id new-id)) (rest lines))]
        (println first-line)
        (doseq [line rest-lines]
          (println line))))))
