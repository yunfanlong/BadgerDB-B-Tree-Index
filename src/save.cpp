/**
 * @author See Contributors.txt for code contributors and overview of BadgerDB.
 * @author Yuinfa Long. Muchan Li
 * @since May 20th, 2023
 * @section LICENSE
 * Copyright (c) 2012 Database Group, Computer Sciences Department, University of Wisconsin-Madison.
 */

#include "btree.h"
#include "filescan.h"
#include "exceptions/bad_index_info_exception.h"
#include "exceptions/bad_opcodes_exception.h"
#include "exceptions/bad_scanrange_exception.h"
#include "exceptions/no_such_key_found_exception.h"
#include "exceptions/scan_not_initialized_exception.h"
#include "exceptions/index_scan_completed_exception.h"
#include "exceptions/file_not_found_exception.h"
#include "exceptions/end_of_file_exception.h"
#include "exceptions/hash_already_present_exception.h"
#include <typeinfo>
#include <iostream>
//#define DEBUG
namespace badgerdb
{

// -----------------------------------------------------------------------------
// BTreeIndex::BTreeIndex -- Constructor
// -----------------------------------------------------------------------------

BTreeIndex::BTreeIndex(const std::string & relationName,
		std::string & outIndexName,
		BufMgr *bufMgrIn,
		const int attrByteOffset,
		const Datatype attrType)
{
  this->scanExecuting = false;
  this->bufMgr = bufMgrIn;
  this->attributeType = attrType;
  this->attrByteOffset = attrByteOffset;

  switch (attrType)
    {
      case INTEGER: 
        this->nodeOccupancy = INTARRAYNONLEAFSIZE;
        this->leafOccupancy = INTARRAYLEAFSIZE;
        break;
      case DOUBLE:
        this->nodeOccupancy = DOUBLEARRAYNONLEAFSIZE;
        this->leafOccupancy = DOUBLEARRAYLEAFSIZE;
        break;
      case STRING:
        this->nodeOccupancy = STRINGARRAYNONLEAFSIZE;
        this->leafOccupancy = STRINGARRAYLEAFSIZE;
        break;
    }

    // Construct the index file name using the relation name and attribute offset
    std::ostringstream idxStr;
    idxStr << relationName << '.' << attrByteOffset;
    std::string indexName = idxStr.str(); // indexName is the name of the index file
    this->outIndexName = indexName; // Return the index file name

    // Check if the index file already exists
    bool indexExists = File::exists(indexName);

    if (indexExists) {
        // Index file exists, open it
        file = new BlobFile(indexName, false); // Open the index file in read/write mode
        this->headerPageNum = file->getFirstPageNo(); // Get the page number of the header page

        // Read the metadata from the header page
        Page* headerPage;
        this->bufMgr->readPage(file, headerPageNum, headerPage);
        IndexMetaInfo* metaInfo = (IndexMetaInfo*)(headerPage);

        // Set the root page number
        this->rootPageNum = metaInfo->rootPageNo;

        // Clean up
        this->bufMgr->unPinPage(file, headerPageNum, false); // Unpin the header page
    } else {
        // Index file does not exist, create it and insert entries from the base relation
        file = new BlobFile(indexName, true); // Create the index file
        this->headerPageNum = file->getFirstPageNo(); // Get the page number of the header page

        // Create the metadata for the index file
        IndexMetaInfo metaInfo;
        strncpy(metaInfo.relationName, relationName.c_str(), sizeof(metaInfo.relationName));
        metaInfo.attrByteOffset = attrByteOffset;
        metaInfo.attrType = attrType;

        Page* rootPage;
        PageId rootPageNo;
        this->bufMgr->allocPage(file, rootPageNo, rootPage);
        this->rootPageNum = rootPageNo;
        this->rootPage = rootPage;

        metaInfo.rootPageNo = rootPageNo;

        Page* emptyPage;
        PageId emptyPageNo;
        this->bufMgr->allocPage(file, emptyPageNo, rootPage);
        this->emptyPageNum = emptyPageNo;
        this->emptyPage = emptyPage;

        Page* nullPage;
        PageId nullPageNo;
        this->bufMgr->allocPage(file, nullPageNo, rootPage);
        this->nullPageNum = nullPageNo;
        this->nullPage = nullPage;


        switch (attrType)
        {
            case INTEGER:
            {
                this->rootNodeInt = reinterpret_cast<NonLeafNode<int>*>(this->rootPage);
                this->rootNodeInt->level = 1;
                this->emptyNodeInt = (LeafNode<int>*) emptyPage;
                this->emptyNodeInt->rightSibPageNo = 0;
                this->rootNodeInt->pageNoArray[0] = emptyPageNum;
                this->nullNodeInt = (NonLeafNode<int>*) nullPage;
                this->nullNodeInt->level = 0;
                break;
            }
            case DOUBLE:
            {
                this->rootNodeDouble = reinterpret_cast<NonLeafNode<double>*>(rootPage);
                this->rootNodeDouble->level = 1;
                this->emptyNodeDouble = (LeafNode<double>*) emptyPage;
                this->emptyNodeDouble->rightSibPageNo = 0;
                this->rootNodeDouble->pageNoArray[0] = emptyPageNum;
                this->nullNodeDouble = (NonLeafNode<double>*) nullPage;
                this->nullNodeDouble->level = 0;
                break;
            }
            case STRING:
            {
                NonLeafNode<std::string>* rootNode = reinterpret_cast<NonLeafNode<std::string>*>(rootPage);
                this->rootNodeStr = reinterpret_cast<NonLeafNode<std::string>*>(rootPage);
                this->rootNodeStr->level = 1;
                this->emptyNodeString = (LeafNode<std::string>*) emptyPage;
                this->emptyNodeString->rightSibPageNo = 0;
                this->rootNodeStr->pageNoArray[0] = emptyPageNum;
                this->nullNodeString = (NonLeafNode<std::string>*) nullPage;
                this->nullNodeString->level = 0;
                break;
            }
            default:
                break;
        }
        this->bufMgr->unPinPage(this->file, emptyPageNum, true);
        this->bufMgr->unPinPage(this->file, nullPageNum, true);
        this->bufMgr->unPinPage(this->file, rootPageNum, true);


        // Write the metadata to the header page
        Page* headerPage;
        bufMgr->allocPage(file, headerPageNum, headerPage);
        memcpy(headerPage, &metaInfo, sizeof(metaInfo));

        this->bufMgr->unPinPage(file, headerPageNum, true); // Unpin the header page

        // Scan the base relation and insert entries into the index
        FileScan fileScan(relationName, bufMgr);
        RecordId rid;
        std::string record;
        while (true) {
            try {
                fileScan.scanNext(rid);
                record = fileScan.getRecord();

                // Extract the key from the record based on the attribute type and offset
                const void* key;
                switch (attrType)
                {
                    case INTEGER:
                        key = (const void*)((const int*)(record.c_str()) + attrByteOffset);
                        // Insert the entry into the index
                        insertEntry(key, rid);
                        break;
                    case DOUBLE:
                        key = (const void*)((const double*)(record.c_str()) + attrByteOffset);
                        insertEntry(key, rid);
                        break;
                    case STRING:
                        key = (const void*)((const char*)(record.c_str()) + attrByteOffset);
                        insertEntry(key, rid);
                        break;
                    default:
                        throw BadIndexInfoException(outIndexName);
                }
            } catch (const EndOfFileException&) {
                // Reached the end of the base relation
                // NonLeafNode<int>* rootNode = (NonLeafNode<int>*) this->rootPage;
                break;
            }
        }
    }
}

// -----------------------------------------------------------------------------
// BTreeIndex::~BTreeIndex -- destructor
// -----------------------------------------------------------------------------
BTreeIndex::~BTreeIndex()
{
    // If there's a scan going on, end it
    try {
        endScan();
    } catch(ScanNotInitializedException e) {
        // Scan already ended or not started, do nothing
    }

    // Write any remaining changes to disk and clear buffer
    bufMgr->flushFile(file);

    // Deallocate the BlobFile object representing the BTree index file
    delete file;
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertEntry
// -----------------------------------------------------------------------------

const void BTreeIndex::insertEntry(const void* key, const RecordId rid)
{
    switch (this->attributeType)
    {
        case INTEGER:
            insertRecursive<int>(this->nullPageNum, this->rootPageNum, *(int*)key, rid);
            break;
        case DOUBLE:
            insertRecursive<double>(this->nullPageNum, this->rootPageNum, *(double*)(key), rid);
            break;
        case STRING:
            char firstTenChar[10];
		        strncpy(firstTenChar, (char *)key, 10);
            insertRecursive<std::string>(this->nullPageNum, this->rootPageNum, firstTenChar, rid);
            break;
        default:
            throw BadIndexInfoException(outIndexName);
    }
}

template <class T>
const void BTreeIndex::insertRecursive(PageId prevPageNo, PageId currentPageNo, T key, RecordId rid)
{
    std::cout << "line 276";
    // Read the data on current page num
    Page* currentPageData;
    this->currentPageNum = currentPageNo;
    this->bufMgr->readPage(this->file, currentPageNo, currentPageData);

    //save the returned page to member var
    this->currentPageData = currentPageData;

    Page* prevPageData;

    this->bufMgr->readPage(this->file, prevPageNo, prevPageData);

    NonLeafNode<T>* prevNode = (NonLeafNode<T>*) prevPageData;

    // Check if the current page is a leaf or non-leaf
    if (prevPageNo != this->nullPageNum && prevNode -> level == 1)
    {
        // at leaf layer
        LeafNode<T>* curLeafNode = reinterpret_cast<LeafNode<T>*>(currentPageData);
        if (curLeafNode->key_count == this->leafOccupancy)
        {
            // Split the leaf node
            PageId newLeafPageNum = splitLeafNode<T>(curLeafNode);
            Page* newLeafPage;
            bufMgr->readPage(this->file, newLeafPageNum, newLeafPage);
			      //create new leaf node / page
            LeafNode<T>* newLeafNode = (LeafNode<T>*)(newLeafPage);

            Page* prevPage;
            this->bufMgr->readPage(this->file, prevPageNo, prevPage);
            T keyToCopyUp;
            copy(keyToCopyUp, newLeafNode->keyArray[0]);
            NonLeafNode<T>* prevNonLeafNode = (NonLeafNode<T>*)(prevPage);
            // copy up the key into NonLeaf Node 
			      insertIntoNonLeafNode<T>(prevNonLeafNode, prevPageNo, keyToCopyUp, currentPageNum, newLeafPageNum);

            // Check if the key should be inserted in the new or old leaf node
            if (compare(key, newLeafNode->keyArray[0])){
              // Insert the key in the new leaf node
              insertIntoLeafNode<T>(newLeafNode, key, rid);
            } else {
              // Insert the key in the old leaf node
              insertIntoLeafNode<T>(curLeafNode, key, rid);
            }

        } else {

            // Leaf node is not full, directly insert the key
          insertIntoLeafNode<T>(curLeafNode, key, rid);
        }
        return;
    }
    else
    {
        
        NonLeafNode<T>* curNonLeafNode = reinterpret_cast<NonLeafNode<T>*>(currentPageData);

		    int orig_count = curNonLeafNode->key_count;

        // Find the next node to follow based on the key
        PageId nextNodePageNum = findPageNoInNonLeaf<T>(curNonLeafNode, key);

        // Recursively insert into the next node
        insertRecursive<T>(this->currentPageNum, nextNodePageNum, key, rid);

        if (curNonLeafNode->key_count == nodeOccupancy)
        {
            // Split the current nonleaf node
            PageId newNonLeafPageNum = splitNonLeafNode<T>(curNonLeafNode);
            // use a new nonleaf page to hold result
            Page* newNonLeafPage;
            this->bufMgr->readPage(this->file, newNonLeafPageNum, newNonLeafPage);
            //create new nonleaf node
            NonLeafNode<T>* newNonLeafNode = (NonLeafNode<T>*)(newNonLeafPage);
            newNonLeafNode->level = curNonLeafNode->level;
            Page* prevPage;
            // if (prevPageNo == this->nullPageNum) {
            // } 
            this->bufMgr->readPage(this->file, prevPageNo, prevPage);
            
            T keyToPushUp;
            copy(keyToPushUp, newNonLeafNode->keyArray[0]);

            copy(newNonLeafNode->keyArray[0] , null_str);

			//TODO: explore index
            for (int i = 0; i < newNonLeafNode->key_count; i++) {
                copy(newNonLeafNode->keyArray[i], newNonLeafNode->keyArray[i + 1]);
				newNonLeafNode->pageNoArray[i] = newNonLeafNode->pageNoArray[i + 1];
            }
            newNonLeafNode->key_count -= 1;
            NonLeafNode<T>* prevNonLeafNode = (NonLeafNode<T>*)(prevPage);
			      insertIntoNonLeafNode<T>(prevNonLeafNode, prevPageNo, keyToPushUp, this->currentPageNum, newNonLeafPageNum);
            this->bufMgr->unPinPage(this->file, currentPageNum, true);
            this->bufMgr->unPinPage(this->file, newNonLeafPageNum, true);
            this->bufMgr->unPinPage(this->file, this->rootPageNum, true);
        }
        else if (curNonLeafNode->key_count > orig_count)
        {
            // Current non-leaf node is not full, the key must have been inserted when processing leaf.
            this->bufMgr->unPinPage(file, this->currentPageNum, true);

        } else {
            this->bufMgr->unPinPage(file, this->currentPageNum, false);
        }

    }// end of else
}

template <class T>
PageId BTreeIndex::splitNonLeafNode(NonLeafNode<T>* nonLeafNode)
{
    // Create a new non-leaf node as the right sibling
    Page* newNonLeafPageData;
    PageId newNonLeafPageNum;
    try {
        this->bufMgr->allocPage(this->file, newNonLeafPageNum, newNonLeafPageData);
    } catch (HashAlreadyPresentException e) {
        
    }
    
    NonLeafNode<T>* newNonLeafNode = (NonLeafNode<T>*)(newNonLeafPageData);

    // Calculate the split point
    int splitIndex = (this->nodeOccupancy) / 2 - 1;
    int iterations = this->nodeOccupancy;

    // Move keys and page numbers to the right sibling
    for (int i = splitIndex; i < iterations; i++) {
        copy(newNonLeafNode->keyArray[i - splitIndex], nonLeafNode->keyArray[i]);
        copy(nonLeafNode->keyArray[i], null_str);
        newNonLeafNode->pageNoArray[i - splitIndex + 1] = nonLeafNode->pageNoArray[i + 1];
        nonLeafNode->pageNoArray[i + 1] = 0;
    }

    //account for the pointer at the beginning
    newNonLeafNode->pageNoArray[0] = nonLeafNode->pageNoArray[splitIndex];
    nonLeafNode->pageNoArray[splitIndex] = 0;

    // Set the key count in the nodes
    nonLeafNode->key_count = splitIndex;
    newNonLeafNode->key_count = iterations - splitIndex;
    
    // Return the page number of the new non-leaf node
    return newNonLeafPageNum;
}

template <class T>
PageId BTreeIndex::splitLeafNode(LeafNode<T>* leafNode)
{
    // Create a new leaf node as the right sibling
    Page* newLeafPageData;
    PageId newLeafPageNum;
    try {
        this->bufMgr->allocPage(file, newLeafPageNum, newLeafPageData);
    } catch (HashAlreadyPresentException e) {

    }
    LeafNode<T>* newLeafNode = (LeafNode<T>*)(newLeafPageData);

    // Calculate the split point
    int splitIndex = (this->leafOccupancy) / 2 - 1;
    int iterations = this->leafOccupancy;

    // Move keys and record IDs to the right sibling
    for (int i = splitIndex; i < iterations; i++) {
        copy(newLeafNode->keyArray[i - splitIndex], leafNode->keyArray[i]);
        //assign the old array positions value to 0 after transferring
        copy(leafNode->keyArray[i] , null_str);
        newLeafNode->ridArray[i - splitIndex] = leafNode->ridArray[i];
    }

    // Set the key count in the nodes
    leafNode->key_count = splitIndex;
    newLeafNode->key_count = iterations - splitIndex;

    // Update the sibling pointers
    newLeafNode->rightSibPageNo = leafNode->rightSibPageNo;
    leafNode->rightSibPageNo = newLeafPageNum;

    // Return the page number of the new leaf node
    return newLeafPageNum;
}

template <class T>
PageId BTreeIndex::findPageNoInNonLeaf(NonLeafNode<T>* nonLeafNode, T key)
{
    if (nonLeafNode->key_count == 0) return nonLeafNode->pageNoArray[0];
    else {
        // Iterate through the keys to find the next node to follow
        for (int i = 0; i <= nonLeafNode->key_count; i++) {
            if (i == nonLeafNode->key_count || !compare(key, nonLeafNode->keyArray[i])) {
                return nonLeafNode->pageNoArray[i];
            }
        }
    }
}


// Insert a new entry into a leaf node
template <class T>
void BTreeIndex::insertIntoLeafNode(LeafNode<T>* leafNode, T& key, RecordId& rid) {
    int insertPos = 0;
    while (insertPos < leafNode->key_count && compare(key, leafNode->keyArray[insertPos]))  {
        insertPos++;
    }
    // Shift existing keys and rids to make space for the new entry
    for (int i = leafNode->key_count; i > insertPos; i--) {
            copy(leafNode->keyArray[i], leafNode->keyArray[i - 1]);
            leafNode->ridArray[i] = leafNode->ridArray[i - 1];
    }
    copy(leafNode->keyArray[insertPos], key);
    leafNode->ridArray[insertPos] = rid;
    leafNode->key_count += 1;
    // }
}

// Insert a new entry into a non-leaf node
template <class T>
void BTreeIndex::insertIntoNonLeafNode(NonLeafNode<T>* nonLeafNode, PageId currPageNo, T key, PageId oldPageNo, PageId newPageNo) {
    int insertPos = 0;
    if (currPageNo == this->nullPageNum) {
        Page* newRootPage;
        PageId newRootPageNum;
        try {
            this->bufMgr->allocPage(this->file, newRootPageNum, newRootPage);
        } catch (HashAlreadyPresentException e) {

        }
        this->rootPageNum = newRootPageNum;
        NonLeafNode<T>* newRootNode = (NonLeafNode<T>*) newRootPage;
        newRootNode->level = 0;
        copy(newRootNode->keyArray[0], key);
        newRootNode->pageNoArray[0] = oldPageNo;
        newRootNode->pageNoArray[1] = newPageNo;
        this->bufMgr->unPinPage(this->file, newRootPageNum, true);
    } else {
        while (insertPos < nonLeafNode->key_count && compare(key, nonLeafNode->keyArray[insertPos]))  {
            insertPos++;
        }
        // Shift existing keys and page numbers to make space for the new entry
        for (int i = nonLeafNode->key_count; i > insertPos; i--) {
            copy(nonLeafNode->keyArray[i], nonLeafNode->keyArray[i - 1]);
            nonLeafNode->pageNoArray[i + 1] = nonLeafNode->pageNoArray[i];
        }
        copy(nonLeafNode->keyArray[insertPos], key);
        nonLeafNode->pageNoArray[insertPos + 1] = newPageNo;
        nonLeafNode->key_count += 1;
    }

}

// Convert to typed keys
template <class T>
T BTreeIndex::convert(const void* key) {
    const T* typedKey = (const T*)(key);
    return *typedKey;
}

template <class T>
bool BTreeIndex::query(Operator lowOp, T lowVal, Operator highOp, T highVal, T currVal) {
    if (lowOp == GT && highOp == LT) {
        return currVal > lowVal && currVal < highVal;
    } else if (lowOp == GT && highOp == LTE) {
        return currVal > lowVal && currVal <= highVal;
    } else if (lowOp == GTE && highOp == LT) {
        return currVal >= lowVal && currVal < highVal;
    } else {
        return currVal >= lowVal && currVal <= highVal;
    }
}

const bool BTreeIndex::compare(int Lvalue, int Rvalue){
    return Lvalue > Rvalue;
}

const bool BTreeIndex::compare(double Lvalue, double Rvalue){
    return Lvalue > Rvalue;
}

const bool BTreeIndex::compare(char Lvalue[], char Rvalue[]){
    return strncmp(Lvalue, Rvalue, 10);
}

const bool BTreeIndex::compare(std::string Lvalue, char Rvalue[]){
    return 1;
}

const void BTreeIndex::copy(int& Lvalue, int Rvalue){
    Lvalue = Rvalue;
}

const void BTreeIndex::copy(double& Lvalue, double Rvalue){
    Lvalue = Rvalue;
}

const void BTreeIndex::copy(char Lvalue[], char Rvalue[]){
    strncpy(Lvalue, Rvalue, 10);
}

const void BTreeIndex::copy(char Lvalue[], std::string Rvalue){
    strncpy(Lvalue, Rvalue.c_str(), 10);
}

const void BTreeIndex::copy(std::string Lvalue, std::string Rvalue){
    Lvalue = Rvalue;
}

const void BTreeIndex::copy(int& Lvalue, std::string Rvalue){
    Lvalue = 0;
}

const void BTreeIndex::copy(double& Lvalue, std::string Rvalue){
    Lvalue = 0.0;
}


const void BTreeIndex::startScan(const void* lowValParm,
				   const Operator lowOpParm,
				   const void* highValParm,
				   const Operator highOpParm) {
	// check operands
	if((lowOpParm != GT && lowOpParm != GTE) || (highOpParm != LT && highOpParm != LTE)) {
		throw BadOpcodesException();
	}
  this->scanExecuting = true;

  this->lowOp = lowOpParm; 
  this->highOp = highOpParm;

	// if lowVal > highVal throw exception
	if(this->attributeType == INTEGER) {
		this->lowValInt = *((int*)lowValParm);
		this->highValInt = *((int*)highValParm);

    if (this->lowValInt > this->highValInt) {
      throw BadScanrangeException();
    }

    PageId curPageNum = this->rootPageNum;
    Page* curPage;
    this->bufMgr->readPage(this->file, curPageNum, curPage);
    NonLeafNode<int>* curPageNode = (NonLeafNode<int>*)curPage;
    PageId nextPageNum;

    while (curPageNode->level != 1){
      nextPageNum = findPageNoInNonLeaf(curPageNode, lowValInt);
      this->bufMgr->readPage(this->file, nextPageNum, curPage);
      curPageNode = (NonLeafNode<int>*) curPage;
    }

    Page* leafPage;
    PageId leafPageNo = findPageNoInNonLeaf(curPageNode, lowValInt);
    this->bufMgr->readPage(this->file, leafPageNo, leafPage);
    LeafNode<int>* leafNode = (LeafNode<int>*)leafPage;

    while (true) {
      for (int i = 0; i < leafNode->key_count; i++){
        if (lowOpParm == GTE && highOpParm == LTE && leafNode->keyArray[i] >= lowValInt && leafNode->keyArray[i] <= highValInt) {
          this->currentPageNum = leafPageNo;
          this->currentPageData = leafPage;
          this->nextEntry = i;

          return;
        } else if (lowOpParm == GT && leafNode->keyArray[i] > lowValInt) {
          this->currentPageNum = leafPageNo;
          this->currentPageData = leafPage;
          this->nextEntry = i;
          return;
        }
      }

        //No such key for the whole loop, throw exception
        if (leafNode->rightSibPageNo == 0){
          this->scanExecuting = false;
          throw NoSuchKeyFoundException();
        }
    }
      
      this->leafPageNum = leafNode->rightSibPageNo;
      this->foundInd = 0;
  }
}