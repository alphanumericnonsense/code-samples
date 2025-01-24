#include <stdio.h>

struct node {
    int label;
    int weight;
    struct node *left;
    struct node *right;
    struct node *parent;
};

void heapify(struct node A[], int i, int heap_size) {
    int smallest;
    int l = 2*i+1;
    int r = 2*i+2;
    if ((l < heap_size) && (A[l].weight < A[i].weight)) {
        smallest = l;
    } else {
        smallest = i;
    }
    if ((r < heap_size) && (A[r].weight < A[smallest].weight)) {
        smallest = r;
    }
    if (smallest != i) {
        struct node tmp = A[i];
        A[i] = A[smallest];
        A[smallest] = tmp;
        heapify(A, smallest, heap_size);
    }
}

void build_heap(struct node A[], int heap_size) {
    for (int i = heap_size / 2 - 1; i >= 0; i--) {
        heapify(A, i, heap_size);
    }
}

struct node extract_min(struct node A[], int heap_size) {
    struct node min = A[0];
    A[0] = A[heap_size - 1]; // lossy...
    heapify(A, 0, heap_size - 1);
    return min;
}

void heap_insert(struct node A[], struct node key, int heap_size) {
    int i = heap_size;
    while ((i > 0) && (A[(i-1)/2].weight > key.weight)) {
        A[i] = A[(i-1)/2];
        i = (i-1)/2;
    }
    A[i] = key;
}

void huffman(struct node A[], int lenA, struct node extracted[]) {
    for (int i = 0; i < lenA; i++) {
        A[i].left = NULL;
        A[i].right = NULL;
        A[i].parent = NULL;
    }

    build_heap(A, lenA);

    struct node internal[lenA - 1];
    int queue_size = lenA;
    for (int i = 0; i < lenA - 1; i++) {
        extracted[2*i] = extract_min(A, queue_size);
        queue_size -= 1;
        extracted[2*i + 1] = extract_min(A, queue_size);
        queue_size -= 1;
        internal[i].weight = extracted[2*i].weight + extracted[2*i+1].weight;
        internal[i].label = -1;
        internal[i].left = &extracted[2*i];
        internal[i].right = &extracted[2*i + 1];
        internal[i].parent = NULL;
        heap_insert(A, internal[i], queue_size);
        queue_size += 1;
    }
    extracted[2*lenA - 2] = A[0];
}

void compute_codes(
    struct node B[],
    struct node *current,
    int code[],
    int array_index,
    int codes[256][256],
    int lengths[]
) {
    if (current->left != NULL) {
        code[array_index] = 0;
        compute_codes(B, current->left, code, array_index + 1, codes, lengths);
    }
    if (current->right != NULL) {
        code[array_index] = 1;
        compute_codes(B, current->right, code, array_index + 1, codes, lengths);
    }
    if (current->label != -1) {
        //printf("%08b: ", current->label);
        lengths[current->label] = array_index;
        for (int i = 0; i < array_index; i++) {
            codes[current->label][i] = code[i];
        }
    }
}
/**********************************************************
 * testing
 **********************************************************/
void test_build_heap() {
    struct node A[15];
    for (int i = 0; i < 15; i++) {
        A[i].weight = (2*i + 3) % 15;
    }
    build_heap(A, 15);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < (1<<i); j++) {
            printf("%d ", A[(1<<i) - 1  + j].weight);
        }
        printf("\n");
    }
    return;
}
/*********************
 * input:
 *
 * 3
 * 5 7
 * 9 11 13 0
 * 2 4 6 8 10 12 14 1
 *
 * output:
 *
 * 0
 * 2 1
 * 4 6 10 3
 * 9 5 11 8 13 12 14 7
 *********************/

/**********************************************************
 * MAIN
 **********************************************************/

// int main() {
//     test_build_heap();
//     return 0;
// }

int main(int argc, char *argv[]) {
    struct node A[256];
    for (int i = 0; i < 256; i++) {
        A[i].label = i;
    }
    struct node B[2*256-1];

    FILE *input_ptr;
    input_ptr = fopen(argv[1], "rb");
    unsigned char byte;
    // printf(argv[1]);
    // printf("\n");

    while (fread(&byte, sizeof(unsigned char), 1, input_ptr) == 1) {
        A[byte].weight += 1;
    }

    fclose(input_ptr);

    int used[256];
    for (int i = 0; i < 256; i++) {
        used[i] = A[i].weight;
    }

    // printf("weights pre-heap:\n");
    // for (int i = 0; i < 8; i++) {
    //     for (int j = 0; j < (1<<i); j++) {
    //         printf("%d ", A[(1<<i) - 1  + j].weight);
    //     }
    //     printf("\n");
    // }

    huffman(A, 256, B);

    int code[256];
    int codes[256][256];
    int lengths[256];
    compute_codes(B, &B[2*256 - 2], code, 0, codes, lengths);

    printf("BIBLE CODES\nUsed bytes\nint, printable ASCII, byte, frequency: code\n\n");

    char x;
    for (int i = 0; i < 256; i++) {
        if (i >= 32 && i <= 126) {
            x = i;
        } else {
            x = 'r';
        }
        if (used[i] > 0) {
            printf("%d %c %08b %d:\t\t", i, x, i, used[i]);
            for (int j = 0; j < lengths[i]; j++) {
                printf("%d", codes[i][j]);
            }
            printf("\n");
        }
    }
    return 0;
}
