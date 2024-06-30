#ifndef ADS_SET_H
#define ADS_SET_H

#include <algorithm>
#include <cassert>
#include <utility>
#include <functional>
#include <iostream>
#include <queue>  //for dump()
#include <vector> //for dump()

#define DEG 30 // d < x < 2d Elements therefore one page has 4 elements

//#define DEBUG
#if defined(DEBUG)
#define TRACE(x) //std::cout << x << "\n";
#define TRACE_DEBUG(x) //std::cout << "DEBUG: " << x << "\n";
#endif

template <typename Key, size_t N = DEG> class ADS_set {
public:
    class Iterator;
    using value_type = Key;
    using key_type = Key;
    using reference = value_type &;
    using const_reference = const value_type &;
    using size_type = size_t;
    using difference_type = std::ptrdiff_t;
    using const_iterator = Iterator;
    using iterator = const_iterator;
    using key_compare = std::less<key_type>;   // B+-Tree
    using key_equal = std::equal_to<key_type>; // Hashing

private:
    // -----------------help structs-----------------------
    struct Leaf;
    struct Internal;
    struct Node; // abstract
    enum class InsState { COPY, PUSH, NONE };

    struct InsertRes {
        key_type copy_key;
        Node *ovf;
        InsState state;
    };

    // -----------------instance variables-----------------
    size_type set_size;
    Node *root;
    int depth;          // only for dump()
    Node *balance_node; //more of a boolean rather than a Node*
    Node *bottom_left;  //iterator

    //---------------private function sigs-----------------

    // --------------private helping functions-------------

    // TODO: find a better way to rewrite this with parent pointer
    InsertRes _insert(const key_type &_key, Node *cur = nullptr,
                      bool redo_search = false) noexcept {
        if (redo_search || count(_key) == 0) {
            // if the key is not found, insert
            if (!cur->is_leaf) {
                // we are not in a leaf, recurse until the end
                Internal *cur_i = static_cast<Internal *>(cur);
                Node *child = cur_i->search_node_internal(_key);
                InsertRes res = _insert(_key, child, true);
                if (res.state != InsState::NONE) {
                    return cur_i->insert_internal(res, root, depth);
                }
                return res;
            }
            set_size++;
            // cur is now a leaf node
            Leaf *cur_l = static_cast<Leaf *>(cur);
            return cur_l->insert_leaf(_key, root, depth);
        }
        // TRACE_DEBUG("SHOULD NOT BE HERE")
        return InsertRes{_key, nullptr, InsState::NONE}; //???
    }

    Node *delete_rebalance(Node *cur, Node *l_neighbor, Node *r_neighbor, Node *l_anc, Node *r_anc,
                           const key_type &key, Node *parent) {

        // clang-format off
        Node *next_child      = nullptr;
        Node *next_l_neighbor = nullptr;
        Node *next_r_neighbor = nullptr;
        Node *next_l_anchor   = nullptr;
        Node *next_r_anchor   = nullptr;
        Node *remove_node     = nullptr;
        uint32_t next_child_index = N * 2 + 1;
        // clang-format on

        if (cur->node_size > N)
            balance_node = nullptr;
        else if (balance_node == nullptr) {
            if (!cur->root_node)
                balance_node = cur;
            else {
                if (cur->node_size <= 1)
                    balance_node = cur;
                else
                    balance_node = nullptr;
            }
        }
        //took me a long ass while to fully get this
        //it keeps track of the topmost node in a certain merge path, meaning that any nodes
        //that go down in this path are also minimally sized

        if (!cur->is_leaf) {
            // we are not in a leaf, recurse
            Internal *i_cur = static_cast<Internal *>(cur);
            next_child_index = i_cur->search_node_st(key);
            next_child = i_cur->child_nodes[next_child_index];
            //find anchors and neighbors
            if (next_child_index == 0) {
                //if were going down completely at the left
                //the left neighbor is going to the rightmost of l_neighbor
                //anchor doesnt change
                Internal *il_neighbor = static_cast<Internal *>(l_neighbor);
                //todo: check correctness
                next_l_neighbor = !l_neighbor
                                      ? nullptr
                                      : il_neighbor->child_nodes[il_neighbor->node_size];
                next_l_anchor = l_anc;
            } else {
                //TODO: check correctness
                next_l_neighbor = i_cur->child_nodes[next_child_index - 1];
                next_l_anchor = cur;
            }
            if (next_child_index == i_cur->node_size) {
                //if were going down completely at the right
                //the right neighbor is going to the leftmost of r_neighbor
                //anchor doesnt change
                Internal *ir_neighbor = static_cast<Internal *>(r_neighbor);
                //TODO: check correctness
                next_r_neighbor = !r_neighbor ? nullptr : ir_neighbor->child_nodes[0];
                next_r_anchor = r_anc;
            } else {
                //TODO: check correctness
                next_r_neighbor = i_cur->child_nodes[next_child_index + 1];
                next_r_anchor = cur;
            }

            //anchors and neighbours are found, go down the recursion path
            remove_node = delete_rebalance(next_child, next_l_neighbor, next_r_neighbor,
                                           next_l_anchor, next_r_anchor, key, cur);
        } else {
            // leaf node
            if (count(key) == 0) {
                balance_node = nullptr;
                return nullptr;
            }
            set_size--;
            assert(cur->is_leaf == true); //we should be in a leaf
            uint32_t i = binary_search(cur->values, cur->node_size, key);
            assert(i!=N*2); //node should exist
            Leaf *l_cur = static_cast<Leaf *>(cur);
            l_cur->remove_leaf_key(i);
            next_child = nullptr;
            remove_node = cur;
        }

        if (remove_node == next_child) {
            //TODO: clear this useless part
            //we need to delete!!!
            assert(cur->is_leaf==false);      //we should be in an internal
            assert(next_child_index < N*2+1); //next_child_index sould not be untouched
            //Internal *i_cur = static_cast<Internal *>(cur);
            //if(remove_node->is_leaf) delete static_cast<Leaf*>(remove_node);
            //else delete static_cast<Internal*>(remove_node);
            //delete remove_node;
            remove_node = nullptr;
            //i_cur->remove_internal_key(next_child_index);
            //note to self: check which ptr we delete here right or left
        }
        Node *to_return;
        if (balance_node == nullptr) {
            to_return = nullptr;
        } else if (cur->root_node) {
            to_return = delete_root(cur);
            //not really useful since delete_root is usually the last call
        } else {
            to_return = balance(cur, l_neighbor, r_neighbor, l_anc, r_anc, parent);
        }
        return to_return;
    }

    [[nodiscard]] Node *delete_root(Node *old_root) {
        if (old_root->is_leaf) //do nothing lol
            return nullptr;
        assert(old_root->node_size==0);

        Internal *i_old_root = static_cast<Internal *>(old_root);
        root = i_old_root->child_nodes[0]; //point to the only child the root will have
        i_old_root->child_nodes[0] = nullptr;
        root->root_node = true;
        balance_node = nullptr; //???
        depth--;
        delete old_root;
        old_root = nullptr;
        return root;
    }

    [[nodiscard]] Node *balance(Node *cur, Node *l_neighbor, Node *r_neighbor, Node *l_anc,
                                Node *r_anc, Node *parent) {
        //choose the best candidate for shifting, one that has more values
        //set the balance_node to that, so we can know if we reached the merge end or not
        //if were merging then use the one that has the closest anchor
        Node *to_return = nullptr; //kinda redundant but i guess copy elision?

        int help = 4;
        assert(help==4);
        if (!l_neighbor) {
            help = r_neighbor->node_size > N;
            to_return = help
                            ? shift(cur, r_neighbor, r_anc, false)
                            : merge(cur, r_neighbor, r_anc, false);
        } else if (!r_neighbor) {
            help = l_neighbor->node_size > N;
            to_return = help
                            ? shift(cur, l_neighbor, l_anc, true)
                            : merge(cur, l_neighbor, l_anc, true);
        } else {
            if (l_anc == parent && r_anc == parent) {
                if (l_neighbor->node_size <= N && r_neighbor->node_size <= N) {
                    help = 1;
                    assert(l_neighbor->node_size == r_neighbor->node_size);
                    to_return = merge(cur, l_neighbor, l_anc, true);
                } else {
                    help = l_neighbor->node_size > r_neighbor->node_size;
                    to_return = help
                                    ? shift(cur, l_neighbor, l_anc, true)
                                    : shift(cur, r_neighbor, r_anc, false);
                }
            } else if (l_anc == parent) {
                help = l_neighbor->node_size > N;
                to_return = help
                                ? shift(cur, l_neighbor, l_anc, true)
                                : merge(cur, l_neighbor, l_anc, true);
            } else if (r_anc == parent) {
                help = r_neighbor->node_size > N;
                to_return = help
                                ? shift(cur, r_neighbor, r_anc, false)
                                : merge(cur, r_neighbor, r_anc, false);
            }
        }

        assert(help!=4);

        return to_return;
    }

    Node *shift(Node *cur, Node *shift_neighb, Node *shift_anc, const bool &l_r) {
        //l_r->true means we got left neighbor, false means we got right
        assert(cur);
        assert(shift_neighb);
        assert(shift_anc);
        Internal *i_cur = static_cast<Internal *>(cur);
        Internal *i_shift_neighb = static_cast<Internal *>(shift_neighb);
        Internal *i_shift_anc = static_cast<Internal *>(shift_anc);
        Leaf *l_cur = static_cast<Leaf *>(cur);
        Leaf *l_shift_neighb = static_cast<Leaf *>(shift_neighb);

        uint32_t anchor_value_index = N * 2 + 1;

        if (!l_r) {
            for (uint32_t i = 0; i < i_shift_anc->node_size; i++) {
                if (i_shift_anc->child_nodes[i] == cur &&
                    i_shift_anc->child_nodes[i + 1] == shift_neighb)
                    anchor_value_index = i;
            }
            if (!cur->is_leaf) {
                //if were gonna shift on an internal, we push down anchors value, and push up the
                //neighbours/this value to anchor
                assert(i_cur->node_size != i_shift_neighb->node_size);
                //well were shifting so logically shouldt be equal
                //we are shifting from right->left
                assert(anchor_value_index<N*2+1);

                i_cur->insert_internal_key(i_shift_anc->values[anchor_value_index],
                                           i_shift_neighb->child_nodes[0]);
                while (i_cur->node_size < i_shift_neighb->node_size - 2) {
                    i_cur->insert_internal_key(i_shift_neighb->values[0],
                                               i_shift_neighb->child_nodes[1]);
                    i_shift_neighb->remove_internal_key(0);
                }
                shift_anc->values[anchor_value_index] = (i_shift_neighb->values[0]);
                i_shift_neighb->remove_internal_key_lp(0);
            } else {
                assert(cur->is_leaf);
                //were in a leaf node, shifting from right->left
                while (l_cur->node_size < l_shift_neighb->node_size) {
                    l_cur->insert_leaf_key(l_shift_neighb->values[0]);
                    l_shift_neighb->remove_leaf_key(0);
                }
                shift_anc->values[anchor_value_index] = shift_neighb->values[0];
                assert(cur->node_size >=N);
                assert(shift_neighb->node_size >=N);
            }

        } else {
            for (uint32_t i = 0; i < i_shift_anc->node_size; i++) {
                if (i_shift_anc->child_nodes[i] == shift_neighb &&
                    i_shift_anc->child_nodes[i + 1] == cur)
                    anchor_value_index = i;
            }
            if (!cur->is_leaf) {
                assert(cur->node_size != shift_neighb->node_size);
                //well were shifting so logically shouldt be equal
                //we are shifting from left->right
                //TODO: unsure about this
                assert(anchor_value_index<N*2+1);
                //TODO: IMPORTANT this condition isnt good and can lead to the neighbor udnerflowing
                i_cur->insert_internal_key_l(i_shift_anc->values[anchor_value_index],
                                             i_shift_neighb->child_nodes[i_shift_neighb->
                                                 node_size]);
                while (i_cur->node_size < i_shift_neighb->node_size - 2) {
                    uint32_t nbsz = i_shift_neighb->node_size;
                    i_cur->insert_internal_key_l(i_shift_neighb->values[nbsz - 1],
                                                 i_shift_neighb->child_nodes[nbsz - 1]);
                    i_shift_neighb->remove_internal_key_lp(nbsz);
                }
                i_shift_anc->values[anchor_value_index] = i_shift_neighb->values[
                    i_shift_neighb->node_size - 1];
                i_shift_neighb->remove_internal_key(i_shift_neighb->node_size - 1);
            } else {
                assert(cur->is_leaf);
                //were in a leaf node, shifting from left->right
                while (l_cur->node_size < l_shift_neighb->node_size) {
                    l_cur->insert_leaf_key(l_shift_neighb->values[l_shift_neighb->node_size - 1]);
                    l_shift_neighb->remove_leaf_key(l_shift_neighb->node_size - 1);
                }
                shift_anc->values[anchor_value_index] = cur->values[0];
                assert(cur->node_size >=N);
                assert(shift_neighb->node_size >=N);
            }
        }

        balance_node = nullptr;
        return nullptr;
    }

    Node *merge(Node *cur, Node *merge_neighb, Node *merge_anc, bool l_r) {
        //l_r==true->neighbor is to the left
        Node *left;
        Node *right;
        if (l_r) {
            left = merge_neighb;
            right = cur;
        } else {
            left = cur;
            right = merge_neighb;
        }
        Internal *i_merge_anc = static_cast<Internal *>(merge_anc);
        Internal *i_left = static_cast<Internal *>(left);
        Internal *i_right = static_cast<Internal *>(right);
        uint32_t anchor_value_index = N * 2 + 1;
        for (uint32_t i = 0; i < i_merge_anc->node_size; i++) {
            //to find slot
            if (i_merge_anc->child_nodes[i] == left &&
                i_merge_anc->child_nodes[i + 1] == right)
                anchor_value_index = i;
        }
        if (!left->is_leaf) {
            assert(anchor_value_index < N*2+1);
            i_left->insert_internal_key(i_merge_anc->values[anchor_value_index],
                                        i_right->child_nodes[0]);
            i_merge_anc->remove_internal_key(anchor_value_index);
            while (right->node_size > 0) {
                i_left->insert_internal_key(i_right->values[0], i_right->child_nodes[1]);
                i_right->remove_internal_key(0);
            }
        } else {
            assert(anchor_value_index < N*2+1);
            Leaf *l_left = static_cast<Leaf *>(left);
            Leaf *l_right = static_cast<Leaf *>(right);
            i_merge_anc->remove_internal_key(anchor_value_index);
            while (right->node_size > 0) {
                l_left->insert_leaf_key(l_right->values[0]);
                l_right->remove_leaf_key(0);
            }
            l_left->next_page = l_right->next_page;
        }
        if (left == balance_node || right == balance_node)
            balance_node = nullptr;
        if (right->is_leaf) delete static_cast<Leaf *>(right);
        else delete static_cast<Internal *>(right);
        right = nullptr;
        return right;
    }

    static uint32_t binary_search(value_type *values, const size_type max, const key_type &_key) {
        if (max == 0)
            return N * 2;
        int left = -1;
        int right = max;                        //???
        int mid = left + ((right - left) >> 1); // avoid overflow

        // TODO: explore possibly branchless, more optimized way to search, maybe unrolling?
        while (right - left > 1) {
            mid = left + ((right - left) >> 1); // avoid , overflow
            if (key_compare{}(_key, values[mid]))
                right = mid;
            else
                left = mid;
        }
        return (left != -1 && key_equal{}(_key, values[left])) ? left : N * 2;
    }

    static int lower_bound(value_type *values, const size_type max, const key_type &_key) {
        int left = 0;
        int right = max - 1;
        int mid = left + ((right - left) >> 1);

        int ans = -1;
        while (left <= right) {
            mid = left + ((right - left) >> 1);
            // the right
            if (!key_compare{}(values[mid], _key))
                right = mid - 1;
            else {
                ans = mid;
                left = mid + 1;
            }
        }
        return ans + 1;
    }

    void delete_tree(Node* cur) {
        if(cur->is_leaf) {
            delete static_cast<Leaf*>(cur);
            return;
        }
        Internal * i_cur = static_cast<Internal*>(cur);
        for(uint32_t i = 0; i < cur->node_size+1; i++) {
            delete_tree(i_cur->child_nodes[i]);
        }
        delete i_cur;

    }

public:
    ADS_set() : set_size{0}, root{new Leaf{true}}, depth{0}, balance_node{nullptr},
                bottom_left{this->root} {
    } // PH1

    ADS_set(std::initializer_list<key_type> ilist) : ADS_set() {
        // PH1
        for (auto &i : ilist)
            _insert(i, root);
    }

    template <typename InputIt> ADS_set(InputIt first, InputIt last) : ADS_set() {
        // PH1
        for (auto &it = first; it != last; it++)
            _insert(*it, root);
    }

    ADS_set(const ADS_set &other) : ADS_set() {
        for (const auto &it : other)
            _insert(it, root, false);
    }

    ~ADS_set() {
        //delete root;
        delete_tree(root);
        root = nullptr;
    }

    ADS_set &operator=(const ADS_set &other) {
        //delete root;
        delete_tree(root);
        root = new Leaf{true};
        set_size = 0;
        depth = 0;
        balance_node = nullptr;
        bottom_left = root;
        //TODO: remember bottom_left
        for (auto &it : other) {
            _insert(it, root, false);
        }
        return *this;
    }

    ADS_set &operator=(std::initializer_list<key_type> ilist) {
        //delete root;
        delete_tree(root);
        root = new Leaf{true};
        set_size = 0;
        depth = 0;
        balance_node = nullptr;
        bottom_left = root;
        //TODO: remember bottom_left
        for (auto &it : ilist) {
            _insert(it, root, false);
        }
        return *this;
    }

    [[nodiscard]] size_type size() const { return set_size; }  // PH1
    [[nodiscard]] bool empty() const { return set_size == 0; } // PH1

    void insert(std::initializer_list<key_type> ilist) {
        // PH1
        for (auto &i : ilist)
            _insert(i, root);
    }

    std::pair<iterator, bool> insert(const key_type &key) {
        //TODO: do whatever this function is
        Iterator it = find(key);
        if (it != end()) return std::make_pair(it, false);
        _insert(key, root, false);
        it = find(key);
        return std::make_pair(it, true);
    }

    template <typename InputIt> void insert(InputIt first, InputIt last) {
        // PH1
        for (auto &it = first; it != last; it++)
            _insert(*it, root);
    }

    void clear() {
        //delete root;
        delete_tree(root);
        root = new Leaf{true};
        set_size = 0;
        depth = 0;
        balance_node = nullptr;
        bottom_left = root;
        //TODO: remember bottom_left
    }

    size_type erase(const key_type &key) {
        size_type test = set_size;
        balance_node = nullptr;
        delete_rebalance(root, nullptr, nullptr, nullptr, nullptr, key, nullptr);
        if (set_size < test) return 1;
        return 0;
    }

    size_type count(const key_type &key) const {
        // PH1
        // Search the whole tree, return 1 when element is found, if not found,
        // return 0

        Node *ptr_old = nullptr;
        Node *ptr = root;
        // IDEA: leaf node will return nullptr if element isnt found,
        // if element found, then Leaf returns a pointer to itsself, stops loop
        while (ptr != nullptr && ptr != ptr_old) {
            ptr_old = ptr;
            if (ptr->is_leaf) {
                Leaf *l_ptr = static_cast<Leaf *>(ptr);
                ptr = l_ptr->search_node_leaf(key);
            } else {
                Internal *i_ptr = static_cast<Internal *>(ptr);
                ptr = i_ptr->search_node_internal(key);
            }
        }
        return ptr == nullptr ? 0 : 1;
    }

    iterator find(const key_type &key) const {

        Node *ptr_old = nullptr;
        Node *ptr = root;
        // IDEA: leaf node will return nullptr if element isnt found,
        // if element found, then Leaf returns a pointer to itsself, stops loop
        while (ptr != nullptr && ptr != ptr_old) {
            ptr_old = ptr;
            if (ptr->is_leaf) {
                Leaf *l_ptr = static_cast<Leaf *>(ptr);
                ptr = l_ptr->search_node_leaf(key);
            } else {
                Internal *i_ptr = static_cast<Internal *>(ptr);
                ptr = i_ptr->search_node_internal(key);
            }
        }
        if (ptr == nullptr) return end();
        uint32_t a = binary_search(ptr->values, ptr->node_size, key);
        assert(ptr->is_leaf);
        Leaf *lptr = static_cast<Leaf *>(ptr);
        return Iterator{lptr, a};
        //TODO: questionable code
    }

    void swap(ADS_set &other) {
        std::swap(root, other.root);
        std::swap(balance_node, other.balance_node);
        std::swap(set_size, other.set_size);
        std::swap(bottom_left, other.bottom_left);
        std::swap(depth, other.depth);
    }

    const_iterator begin() const {
        if (set_size == 0) return end();
        /*Node *ptr = root;
        while (!ptr->is_leaf)
            ptr = static_cast<Internal *>(ptr)->child_nodes[0];
        assert(ptr->is_leaf);
        Leaf *lptr = static_cast<Leaf *>(ptr);
        return Iterator{lptr, 0};*/
        return Iterator{static_cast<Leaf *>(bottom_left), 0};
    }

    const_iterator end() const {
        return Iterator{nullptr, 0};
    }

    void dump(std::ostream &o = std::cerr, int dp = -1) const {
        if (dp == -1)
            dp = depth;
        std::queue<Node *> bfs;
        std::vector<int> v(4, 0);
        v[0] = 1;
        bfs.push(root);
        int cnt = 0;
        while (!bfs.empty()) {
            Node *pop = bfs.front();
            bfs.pop();
            if (!pop->is_leaf) {
                Internal *pop_i = static_cast<Internal *>(pop);
                for (size_type i = 0; i < pop_i->node_size + 1; i++) {
                    bfs.push(pop_i->child_nodes[i]);
                    v[cnt + 1]++;
                }
            }

            for (int i = 0; i < dp; i++)
                o << "\t";
            o << "[" << pop->values[0];
            for (size_type i = 1; i < pop->node_size; i++)
                o << ", " << pop->values[i];
            // if(pop->is_leaf() && static_cast<Leaf*>(pop)->next_page == bfs.front())
            // o << "]" << pop->node_size << "->";
            // else
            o << "]" << pop->node_size << (pop->root_node ? "r " : "  ");

            v[cnt]--;

            if (v[cnt] == 0) {
                o << "\n\n";
                dp--;
                cnt++;
            }
        }
        o << "\n\n depth: " << depth << "\n";
    }

    friend bool operator==(const ADS_set &lhs, const ADS_set &rhs) {
        if (lhs.set_size != rhs.set_size) return false;
        for (auto &it : lhs) {
            if (rhs.count(it) == 0)
                return false;
        }
        return true;
    }

    friend bool operator!=(const ADS_set &lhs, const ADS_set &rhs) {
        return !(rhs == lhs);
    }

private:
    struct Node {
        key_type values[N * 2];
        uint32_t node_size;
        bool root_node; //TODO: remove this
        bool is_leaf;

        explicit Node(const bool root = false)
            : values{}, node_size{0}, root_node{root}, is_leaf{false} {
        }

        ~Node() = default;

    };

    // TODO: measure if final has any performance benefits in these structs
    struct Leaf : Node {
        Leaf *next_page;

        explicit Leaf(const bool root = false, Leaf *next_page = nullptr)
            : Node(root), next_page{next_page} {
            this->is_leaf = true;
        }

        Node *search_node_leaf(const key_type &search_key) noexcept {
            // Run a binary search in the nodes aa
            // if not found, nullptr, otherwise return pointer to itsself
            uint32_t res = binary_search(this->values, this->node_size, search_key);
            if (res >= N * 2)
                return nullptr;
            return this;
        }

        InsertRes insert_leaf(const key_type &_key, Node *&root,
                              int &depth) noexcept {
            if (this->node_size < N * 2) {
                this->insert_leaf_key(_key);
                // TRACE_DEBUG(_key << " is added normally");
                return InsertRes{(_key), nullptr, InsState::NONE};
            }
            if (this->node_size >= N * 2) {
                // deal with overflow
                // my method here is linearly copying elements
                // has room for optimization

                // TODO: here i do the linked list thing, check if right or not
                Leaf *overflow_node = new Leaf(false);
                //TRACE_DEBUG("NEW CALLED")
                overflow_node->next_page = this->next_page;
                this->next_page = overflow_node;

                bool lr = !key_compare{}(_key, this->values[N - 1]) ? true : false;
                // HERE: choose in which node to put the goofy ahh _key, when lr is true

                for (uint32_t i = lr ? N : N - 1; i < this->node_size;) {
                    overflow_node->insert_leaf_key(this->values[i]);
                    this->remove_leaf_key(i);
                }
                // insert the remaining _key
                if (lr) overflow_node->insert_leaf_key(_key);
                else this->insert_leaf_key(_key);

                if (this->root_node) {
                    // cur is root, create new parent node
                    depth++;
                    Internal *ovf_root = new Internal(true);
                    this->root_node = false;
                    root = ovf_root;
                    ovf_root->insert_internal_root(overflow_node->values[0], this,
                                                   overflow_node);
                    return InsertRes{_key, overflow_node, InsState::NONE};
                }
                // current best idea: struct that returns what we should do with the
                // key, node
                return InsertRes{(overflow_node->values[0]), overflow_node,
                                 InsState::COPY};
            }
            // TRACE_DEBUG("THIS CONTROL PATH SHOUL LOGICALLY NOT BE REACHED")
            return InsertRes{_key, nullptr, InsState::NONE};
        }

        void insert_leaf_key(const key_type &key) noexcept {
            // NOTE: this function assumes bounds are already checked
            uint32_t i = static_cast<uint32_t>(lower_bound(this->values, this->node_size, key));
            ++this->node_size;
            for (size_type j = this->node_size - 1; j > i; j--)
                this->values[j] = (this->values[j - 1]);
            this->values[i] = (key);
        }

        void remove_leaf_key(size_type index) noexcept {
            if (this->node_size > 0) this->node_size--;
            for (uint32_t i = index; i < this->node_size; i++) {
                this->values[i] = (this->values[i + 1]);
            }
        }

        ~Leaf() = default;
    };

    // final keyword shows that Internal/Leaf nodes wont have anymore children
    struct Internal : Node {
        //TODO:final
        Node *child_nodes[N * 2 + 1]; // array of Pointers to child nodes

        explicit Internal(const bool root = false) : Node(root), child_nodes{} {
            this->is_leaf = false;
        }

        Node *search_node_internal(const key_type &search_key) noexcept {
            // Run a binary search in the nodes
            // should find the lower bound of this search key
            // find the index of that, then the pointer should simply be +1 from that index
            int ans = lower_bound(this->values, this->node_size, search_key);
            // in case key is equal
            if (static_cast<size_type>(ans) < this->node_size &&
                key_equal{}(this->values[ans], search_key))
                return child_nodes[ans + 1];
            return child_nodes[ans];
        }

        int search_node_st(const key_type &search_key) noexcept {
            //TODO:REMOVE THE OTHER SEARCH AND USE THIS ONLY BECAUSE ITS STUPID
            int ans = lower_bound(this->values, this->node_size, search_key);
            // in case key is equal
            if (static_cast<uint32_t>(ans) < this->node_size &&
                key_equal{}(this->values[ans], search_key))
                return ans + 1;
            return ans;
        }

        InsertRes insert_internal(InsertRes res, Node *&root, int &depth) noexcept {
            InsertRes _res{res.copy_key, nullptr, InsState::NONE};
            // idea here is to insert middle key to the parent
            if (this->node_size < N * 2) {
                this->insert_internal_key(res.copy_key, res.ovf);
                // TRACE_DEBUG(res.copy_key << " is inserted internal normally")
                //return InsertRes{(res.copy_key), nullptr, InsState::NONE};
                return _res;
            }
            if (this->node_size >= N * 2) {
                bool lr =
                    !key_compare{}(res.copy_key, this->values[N - 1]) ? true : false;
                // HERE: choose in which node to put the goofy ahh _key, when lr is true
                // go right overflow_node->insert_leaf_key(_key);

                Internal *overflow_int = new Internal(false);
                //TRACE_DEBUG("NEW CALLED")
                for (uint32_t i = lr ? N : N - 1; i < this->node_size;) {
                    overflow_int->insert_internal_key(this->values[i],
                                                      this->child_nodes[i + 1]);
                    this->remove_internal_key(i);
                }

                // insert the remaining _key
                if (lr)
                    overflow_int->insert_internal_key(res.copy_key, res.ovf);
                else
                    this->insert_internal_key(res.copy_key, res.ovf);
                // now we should get ovf[0] remove it without pointer, and push it up
                // with a pointer to overflow_int
                key_type _k = overflow_int->remove_internal_key_p(0);
                for (uint32_t i = 0; i < overflow_int->node_size + 1; i++)
                    overflow_int->child_nodes[i] = overflow_int->child_nodes[i + 1];
                if (this->root_node) {
                    depth++;
                    Internal *ovf_root = new Internal(true);
                    this->root_node = false;
                    ovf_root->insert_internal_root(_k, this, overflow_int);
                    root = ovf_root;
                    _res.copy_key = _k;
                    _res.ovf = nullptr;
                    _res.state = InsState::NONE;
                    return _res;
                    //return InsertRes{(_k), nullptr, InsState::NONE};
                }
                _res.copy_key = _k;
                _res.ovf = overflow_int;
                _res.state = InsState::PUSH;
                return _res;
                //return InsertRes{_k, overflow_int, InsState::PUSH};
            }
            //TRACE_DEBUG("THIS CONTROL PATH SHOUL LOGICALLY NOT BE REACHED")
            return InsertRes{(res.copy_key), nullptr, InsState::NONE};
        }

        //this funciton probably not needed
        void insert_internal_key_l(const key_type &key, Node *ptr) noexcept {
            ++this->node_size;
            for (uint32_t i = this->node_size; i > 0; i--) {
                this->values[i] = this->values[i - 1];
            }
            for (size_type i = this->node_size + 1; i > 0; i--) {
                this->child_nodes[i] = this->child_nodes[i - 1];
            }
            this->values[0] = key;
            this->child_nodes[0] = ptr;
        }

        void insert_internal_key(const key_type &_key, Node *ptr) noexcept {
            // NOTE: this function assumes bounds are already checked
            uint32_t i = static_cast<uint32_t>(lower_bound(this->values, this->node_size, _key));
            ++this->node_size;
            for (uint32_t j = this->node_size - 1; j > i; j--) {
                this->values[j] = (this->values[j - 1]);
                // TODO: check this, it works but im not sure about it
                this->child_nodes[j + 1] = this->child_nodes[j];
                // avoid off-by-1 error here i dont remember
            }
            this->values[i] = _key;
            this->child_nodes[i + 1] = ptr;
        }

        void insert_internal_root(const key_type &_key, Node *before_ptr,
                                  Node *after_ptr) noexcept {
            this->node_size++;
            this->values[0] = _key;
            this->child_nodes[0] = before_ptr;
            this->child_nodes[1] = after_ptr;
        }

        void remove_internal_key(uint32_t i) noexcept {
            // remove key with right pointer
            if (this->node_size > 0) this->node_size--;
            for (size_type j = i; j < this->node_size; j++) {
                // TODO: check for OFFBY1 error here its prob wrong
                this->values[j] = (this->values[j + 1]);
                this->child_nodes[j + 1] = this->child_nodes[j + 2];
            }

        }

        void remove_internal_key_lp(uint32_t i) noexcept {
            // remove key with left pointer
            if (this->node_size > 0) {
                for (size_type j = i; j < this->node_size - 1; j++) {
                    // TODO: check for OFFBY1 error here its prob wrong
                    this->values[j] = this->values[j + 1];
                }
                for (size_type j = i; j < this->node_size; j++) {
                    // TODO: check for OFFBY1 error here its prob wrong
                    this->child_nodes[j] = this->child_nodes[j + 1];
                }
                this->node_size--;
            }
        }

        key_type remove_internal_key_p(uint32_t i) noexcept {
            // remove a key wihtout a pointer (when we push middle value)
            // redundant code, can just pass a bool in the param, currently just need
            // to get it working
            if (this->node_size > 0) this->node_size--;

            key_type copy = (this->values[i]);
            for (size_type j = i; j < this->node_size; j++) {
                // TODO: check for OFFBY1 error here its prob wrong as
                this->values[j] = (this->values[j + 1]);
            }

            return copy;
        }

        /*~Internal() override {
            if (this->node_size > 0) {
                for (uint32_t i = 0; i < this->node_size + 1; i++) {
                    if (child_nodes[i]->is_leaf)
                        delete static_cast<Leaf *>(child_nodes[i]);
                    else
                        delete static_cast<Internal *>(child_nodes[i]);
                }
            }
        }*/ //To remove vptr
        ~Internal() = default;
    };
};

template <typename Key, size_t N>
class ADS_set<Key, N>::Iterator {
    Leaf *node_ptr;
    uint32_t index;

public:
    using value_type = Key;
    using difference_type = std::ptrdiff_t;
    using reference = const value_type &;
    using pointer = const value_type *;
    using iterator_category = std::forward_iterator_tag;

    explicit Iterator(Leaf *ptr = nullptr,
                      const uint32_t index = 0) : node_ptr{ptr}, index{index} {
    }

    reference operator*() const {
        return node_ptr->values[index];
    }

    pointer operator->() const {
        return node_ptr->values + index;
    }

    Iterator &operator++() {
        if (node_ptr == nullptr) return *this;
        if (index < node_ptr->node_size - 1) {
            ++index;
            return *this;
        }
        node_ptr = node_ptr->next_page;
        index = 0;
        return *this;
        //assert(node_ptr==nullptr); //in this case were at the end
        //index = 0;
        //return *this;
    }

    Iterator operator++(int) {
        Iterator copy{node_ptr, index};
        ++(*this);
        return copy;
    }

    friend bool operator==(const Iterator &lhs, const Iterator &rhs) {
        if (lhs.node_ptr == rhs.node_ptr && lhs.index == rhs.index)
            return true;
        return false;
    }

    friend bool operator!=(const Iterator &lhs, const Iterator &rhs) {
        return !(lhs == rhs);
    }
};


template <typename Key, size_t N>
void swap(ADS_set<Key, N> &lhs, ADS_set<Key, N> &rhs) { lhs.swap(rhs); }


//template <typename Key, size_t N> ADS_set<Key, N>::Node::~Node() {}

#endif // ADS_SET_H