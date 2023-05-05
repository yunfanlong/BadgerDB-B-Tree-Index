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
#include "exceptions/page_not_pinned_exception.h"
#include "exceptions/page_pinned_exception.h"
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

    //fprintf(stdout, "line 43");


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
        this->file = new BlobFile(indexName, false); // Open the index file in read/write mode
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
        this->file = new BlobFile(indexName, true); // Create the index file
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
        this->bufMgr->allocPage(file, emptyPageNo, emptyPage);
        this->emptyPageNum = emptyPageNo;
        this->emptyPage = emptyPage;

        Page* nullPage;
        PageId nullPageNo;
        this->bufMgr->allocPage(file, nullPageNo, nullPage);
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

                break;
            }
            case DOUBLE:
            {
                this->rootNodeDouble = reinterpret_cast<NonLeafNode<double>*>(rootPage);
                this->rootNodeDouble->level = 1;
                this->emptyNodeDouble = (LeafNode<double>*) emptyPage;
                this->emptyNodeDouble->rightSibPageNo = 0;
                this->rootNodeDouble->pageNoArray[0] = emptyPageNum;

                break;
            }
            case STRING:
            {
                this->rootNodeString = reinterpret_cast<NonLeafNode<std::string>*>(rootPage);
                this->rootNodeString  = reinterpret_cast<NonLeafNode<std::string>*>(rootPage);
                this->rootNodeString ->level = 1;
                this->emptyNodeString = (LeafNode<std::string>*) emptyPage;
                this->emptyNodeString->rightSibPageNo = 0;
                this->rootNodeString ->pageNoArray[0] = emptyPageNum;

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

        int count = 0;
        while (true) {
            try {
                fileScan.scanNext(rid);
                record = fileScan.getRecord();

                // Extract the key from the record based on the attribute type and offset
                const void* key;
                key = (const void*)(record.c_str() + attrByteOffset);
                // Insert the entry into the index
                insertEntry(key, rid);
                count++;

            } catch (const EndOfFileException&) {
                // Reached the end of the base relation
                break;
            }
      }
      this->scanExecuting = false;
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
    this->bufMgr->flushFile(this->file);

    // Deallocate the BlobFile object representing the BTree index file
    delete this->file;
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertEntry
// -----------------------------------------------------------------------------

const void BTreeIndex::insertEntry(const void* key, const RecordId rid)
{
    //std::cout <<"Initial Key:  " << key << "/  ";
    if (this->attributeType == INTEGER){
        insertRecursive<int>(this->nullPageNum, this->rootPageNum, *(int*)key, rid);
        return;
    } else if (this->attributeType == DOUBLE){
        insertRecursive<double>(this->nullPageNum, this->rootPageNum, *(double*)(key), rid);
        return;
    } else if (this->attributeType == STRING){
      std::string firstTenChar = std::string((char *)key).substr(0, 10);
      // std::cout << "First Ten Char:  " << firstTenChar << "/  ";
      insertRecursive<std::string>(this->nullPageNum, this->rootPageNum, firstTenChar, rid);
      return;
    } else {
      throw BadIndexInfoException(this->outIndexName);
    }
}

template <class T>
const void BTreeIndex::insertRecursive(PageId prevPageNo, PageId currentPageNo, T key, RecordId rid)
{
    // Read the current page into memory
    Page* currentPageData;
    // Page* prevPageData;
    this->bufMgr->readPage(this->file, currentPageNo, currentPageData);

    Page* prevPageData;

    this->bufMgr->readPage(this->file, prevPageNo, prevPageData);

    NonLeafNode<T>* prevNode = (NonLeafNode<T>*) prevPageData;

    // Check if the current page is a leaf or non-leaf
    if (prevPageNo != this->nullPageNum && prevNode -> level == 1)
    {
        // at leaf layer
        LeafNode<T>* curLeafNode = reinterpret_cast<LeafNode<T>*>(currentPageData);
        if (curLeafNode->key_count == leafOccupancy)
        {
            // Split the leaf node
            PageId newLeafPageNum = splitLeafNode<T>(curLeafNode);
            Page* newLeafPage;
            this->bufMgr->readPage(this->file, newLeafPageNum, newLeafPage);
			      //create new leaf node / page
            LeafNode<T>* newLeafNode = (LeafNode<T>*)(newLeafPage);

            //Retrieve Key to Copy Up
            T keyToCopyUp;
            copy(keyToCopyUp, newLeafNode->keyArray[0]);
            //Retrieve Prev Page Node for Copy Up
            NonLeafNode<T>* prevNonLeafNode = (NonLeafNode<T>*)(prevPageData);
            //Insert the Key into Previous Non Leaf Node
			      insertIntoNonLeafNode<T>(prevNonLeafNode, prevPageNo, keyToCopyUp, currentPageNo, newLeafPageNum);
            //Previous Page Modified, Unpin and mark dirty.


            this->bufMgr->unPinPage(this->file, prevPageNo, true);

            // Check if the key should be inserted in the new or old leaf node
            if (compare(key, newLeafNode->keyArray[0]))
            {
                // Insert the key in the new leaf node
                insertIntoLeafNode<T>(newLeafNode, key, rid);
            }
            else
            {
                // Insert the key in the old leaf node
                insertIntoLeafNode<T>(curLeafNode, key, rid);
            }

            this->bufMgr->unPinPage(this->file, newLeafPageNum, true);
        }
        else
        {
            // Leaf node is not full, directly insert the key
            insertIntoLeafNode<T>(curLeafNode, key, rid);
            this->bufMgr->unPinPage(this->file, prevPageNo, false);
        }
        this->bufMgr->unPinPage(this->file, currentPageNo, true);
        return;
    }
    else
    {
        NonLeafNode<T>* curNonLeafNode = reinterpret_cast<NonLeafNode<T>*>(currentPageData);

		    int orig_count = curNonLeafNode->key_count;

        // Find the next node to follow based on the key
        PageId nextNodePageNum = findPageNoInNonLeaf<T>(curNonLeafNode, key);

        // Recursively insert into the next node
        insertRecursive<T>(currentPageNo, nextNodePageNum, key, rid);

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
            T keyToPushUp;
            copy(keyToPushUp, newNonLeafNode->keyArray[0]);
            copy(newNonLeafNode->keyArray[0] , 0);



            for (int i = 0; i < newNonLeafNode->key_count; i++) {
                copy(newNonLeafNode->keyArray[i], newNonLeafNode->keyArray[i + 1]);
				        newNonLeafNode->pageNoArray[i] = newNonLeafNode->pageNoArray[i + 1];
            }
            newNonLeafNode->key_count -= 1;
            NonLeafNode<T>* prevNonLeafNode = (NonLeafNode<T>*)(prevPageData);
			      insertIntoNonLeafNode<T>(prevNonLeafNode, prevPageNo, keyToPushUp, currentPageNo, newNonLeafPageNum);
            this->bufMgr->unPinPage(this->file, currentPageNo, true);
            this->bufMgr->unPinPage(this->file, newNonLeafPageNum, true);
            this->bufMgr->unPinPage(this->file, prevPageNo, true);
        }
        else if (curNonLeafNode->key_count > orig_count)
        {
            // Current non-leaf node is not full, directly insert the key
            this->bufMgr->unPinPage(file, currentPageNo, true);

        } else {
            this->bufMgr->unPinPage(file, currentPageNo, false);
        }
        //unpin prev page as normally would
        this->bufMgr->unPinPage(this->file, prevPageNo, false);

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
    int splitIndex = nonLeafNode->key_count / 2;
    int iterations = nonLeafNode->key_count;

    // Move keys and page numbers to the right sibling
    for (int i = splitIndex; i < iterations; i++) {
        copy(newNonLeafNode->keyArray[i - splitIndex], nonLeafNode->keyArray[i]);
        copy(nonLeafNode->keyArray[i], 0);
        newNonLeafNode->pageNoArray[i - splitIndex + 1] = nonLeafNode->pageNoArray[i + 1];
        nonLeafNode->pageNoArray[i + 1] = 0;
    }

    //account for the pointer at the beginning
    newNonLeafNode->pageNoArray[0] = nonLeafNode->pageNoArray[splitIndex];
    nonLeafNode->pageNoArray[splitIndex] = 0;

    // Set the key count in the nodes
    nonLeafNode->key_count = splitIndex;
    newNonLeafNode->key_count = iterations - splitIndex;

    this->bufMgr->unPinPage(this->file, newNonLeafPageNum, true);
    
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
    int splitIndex = leafNode->key_count / 2;
    int iterations = leafNode->key_count;

    // Move keys and record IDs to the right sibling
    for (int i = splitIndex; i < iterations; i++) {
        copy(newLeafNode->keyArray[i - splitIndex], leafNode->keyArray[i]);
        copy(leafNode->keyArray[i] , 0);
        newLeafNode->ridArray[i - splitIndex] = leafNode->ridArray[i];
    }

    // Set the key count in the nodes
    leafNode->key_count = splitIndex;
    newLeafNode->key_count = iterations - splitIndex;

    // Update the sibling pointers
    newLeafNode->rightSibPageNo = leafNode->rightSibPageNo;
    leafNode->rightSibPageNo = newLeafPageNum;

    this->bufMgr->unPinPage(this->file, newLeafPageNum, true);

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
        PageId newRootNum;
        try {
            this->bufMgr->allocPage(this->file, newRootNum, newRootPage);
        } catch (HashAlreadyPresentException e) {

        }
        this->rootPageNum = newRootNum;
        NonLeafNode<T>* newRootNode = (NonLeafNode<T>*) newRootPage;
        newRootNode->level = 0;
        copy(newRootNode->keyArray[0], key);
        newRootNode->pageNoArray[0] = oldPageNo;
        newRootNode->pageNoArray[1] = newPageNo;
        this->bufMgr->unPinPage(this->file, newRootNum, true);
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

const bool BTreeIndex::compare(int Lvalue, int Rvalue){
    return Lvalue > Rvalue;
}

const bool BTreeIndex::compare(double Lvalue, double Rvalue){
    return Lvalue > Rvalue;
}


const bool BTreeIndex::compare(std::string& Lvalue, char Rvalue[]){
    return Lvalue > std::string(Rvalue).substr(0, 10);
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

const void BTreeIndex::copy(std::string& Lvalue, char Rvalue[]){
    Lvalue = std::string(Rvalue).substr(0, 10);
}

const void BTreeIndex::copy(char Lvalue[], std::string Rvalue){
    strncpy(Lvalue, Rvalue.c_str(), 10);
}

const void BTreeIndex::copy(std::string Lvalue, std::string Rvalue){
    Lvalue = Rvalue;
}

const void BTreeIndex::copy(char Lvalue[], int Rvalue){
    strncpy(Lvalue, null_str, 10);
}



// -----------------------------------------------------------------------------
// BTreeIndex::startScan
// -----------------------------------------------------------------------------

const void BTreeIndex::startScan(const void* lowValParm,
				   const Operator lowOpParm,
				   const void* highValParm,
				   const Operator highOpParm) {


	if(lowOpParm != GT && lowOpParm != GTE) {
		throw BadOpcodesException();
	}

  if(highOpParm != LT && highOpParm != LTE) {
    throw BadOpcodesException();
  }

  this->lowOp = lowOpParm; 
  this->highOp = highOpParm;

	// if lowVal > highVal throw exception
	if(this->attributeType == INTEGER) {
		this->lowValInt = *((int*)lowValParm);
		this->highValInt = *((int*)highValParm);

    if (this->lowValInt > this->highValInt) {
      throw BadScanrangeException();
    }

    this->scanExecuting = true;

    PageId curPageNum = this->rootPageNum;
    Page* curPage;
    this->bufMgr->readPage(this->file, curPageNum, curPage);
    this->curNodeInt = (NonLeafNode<int>*)curPage;

    // traverse until one level above leaf
    while (this->curNodeInt->level != 1){
      this->bufMgr->unPinPage(this->file, curPageNum, false);
      curPageNum = findPageNoInNonLeaf<int>(this->curNodeInt, lowValInt);
      this->bufMgr->readPage(this->file, curPageNum, curPage);
      this->curNodeInt = (NonLeafNode<int>*) curPage;
    }

    this->bufMgr->unPinPage(this->file, curPageNum, false);

    // find the leaf node that is at least greater than lowVal
    Page* leafPage;
    PageId leafPageNo = findPageNoInNonLeaf<int>(this->curNodeInt, lowValInt);
    this->bufMgr->readPage(this->file, leafPageNo, leafPage);
    LeafNode<int>* leafNode = (LeafNode<int>*)leafPage;
    this->bufMgr->unPinPage(this->file, leafPageNo, false);

    // move to the right neighbor if all keys of the current node are
    // not good
    for (int i = 0; i < leafNode->key_count; i++){
      if (lowOpParm == GTE && leafNode->keyArray[i] >= lowValInt) {
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
      if (leafNode->rightSibPageNo != 0){
        this->currentPageNum = leafNode->rightSibPageNo;
        this->bufMgr->readPage(this->file, this->currentPageNum, this->currentPageData);
        this->nextEntry = 0;
        this->bufMgr->unPinPage(this->file, this->currentPageNum, false);
        return;
      }
      
      this->scanExecuting = false;
      throw NoSuchKeyFoundException();
  } else if(this->attributeType == DOUBLE) { // same as detailed in the int case
      this->lowValDouble = *((double*)lowValParm);
      this->highValDouble = *((double*)highValParm);

      if (this->lowValDouble > this->highValDouble) {
        throw BadScanrangeException();
      }

      this->scanExecuting = true;

      PageId curPageNum = this->rootPageNum;
      Page* curPage;
      this->bufMgr->readPage(this->file, curPageNum, curPage);
      this->curNodeDouble = (NonLeafNode<double>*)curPage;

      while (this->curNodeDouble->level != 1){
        this->bufMgr->unPinPage(this->file, curPageNum, false);
        curPageNum = findPageNoInNonLeaf<double>(this->curNodeDouble, lowValDouble);
        this->bufMgr->readPage(this->file, curPageNum, curPage);
        this->curNodeDouble = (NonLeafNode<double>*) curPage;
      }

      this->bufMgr->unPinPage(this->file, curPageNum, false);

      Page* leafPage;
      PageId leafPageNo = findPageNoInNonLeaf<double>(this->curNodeDouble, lowValDouble);
      this->bufMgr->readPage(this->file, leafPageNo, leafPage);
      LeafNode<double>* leafNode = (LeafNode<double>*)leafPage;
      this->bufMgr->unPinPage(this->file, leafPageNo, false);

      for (int i = 0; i < leafNode->key_count; i++){
        if (lowOpParm == GTE && leafNode->keyArray[i] >= lowValDouble) {
          this->currentPageNum = leafPageNo;
          this->currentPageData = leafPage;
          this->nextEntry = i;
          return;
        } else if (lowOpParm == GT && leafNode->keyArray[i] > lowValDouble) {
          this->currentPageNum = leafPageNo;
          this->currentPageData = leafPage;
          this->nextEntry = i;
          return;
        }
      }

      //No such key for the whole loop, throw exception
      if (leafNode->rightSibPageNo != 0){
        this->currentPageNum = leafNode->rightSibPageNo;
        this->bufMgr->readPage(this->file, this->currentPageNum, this->currentPageData);
        this->nextEntry = 0;
        this->bufMgr->unPinPage(this->file, this->currentPageNum, false);
        return;
      }
      
      this->scanExecuting = false;
      throw NoSuchKeyFoundException();
  } else if(this->attributeType == STRING) { // same as detailed in the int case
      this->lowValString = std::string((char*)lowValParm).substr(0, 10);
      this->highValString = std::string((char*)highValParm).substr(0, 10);;

      if (this->lowValString > this->highValString) {
        throw BadScanrangeException();
      }

      this->scanExecuting = true;

      PageId curPageNum = this->rootPageNum;
      Page* curPage;
      this->bufMgr->readPage(this->file, curPageNum, curPage);
      this->curNodeString = (NonLeafNode<std::string>*)curPage;

      while (this->curNodeString->level != 1){
        this->bufMgr->unPinPage(this->file, curPageNum, false);
        curPageNum = findPageNoInNonLeaf<std::string>(this->curNodeString, lowValString);
        this->bufMgr->readPage(this->file, curPageNum, curPage);
        this->curNodeString = (NonLeafNode<std::string>*) curPage;
      }

      this->bufMgr->unPinPage(this->file, curPageNum, false);

      Page* leafPage;
      PageId leafPageNo = findPageNoInNonLeaf<std::string>(this->curNodeString, lowValString);
      this->bufMgr->readPage(this->file, leafPageNo, leafPage);
      LeafNode<std::string>* leafNode = (LeafNode<std::string>*)leafPage;
      this->bufMgr->unPinPage(this->file, leafPageNo, false);

      for (int i = 0; i < leafNode->key_count; i++){
        if (lowOpParm == GTE && (std::string(leafNode->keyArray[i]).substr(0, 10) >= lowValString)) {
          this->currentPageNum = leafPageNo;
          this->currentPageData = leafPage;
          this->nextEntry = i;
          return;
        } else if (lowOpParm == GT && (std::string(leafNode->keyArray[i]).substr(0, 10) > lowValString)) {
          this->currentPageNum = leafPageNo;
          this->currentPageData = leafPage;
          this->nextEntry = i;
          return;
        }
      }

      //No such key for the whole loop, throw exception
      if (leafNode->rightSibPageNo != 0){
        this->currentPageNum = leafNode->rightSibPageNo;
        this->bufMgr->readPage(this->file, this->currentPageNum, this->currentPageData);
        this->nextEntry = 0;
        this->bufMgr->unPinPage(this->file, this->currentPageNum, false);
        return;
      }
      
      this->scanExecuting = false;
      // this->bufMgr->unPinAllPages(this->file);
      throw NoSuchKeyFoundException();
  }
} 






// -----------------------------------------------------------------------------
// BTreeIndex::scanNext
// -----------------------------------------------------------------------------

const void BTreeIndex::scanNext(RecordId& outRid)
{
  if(!scanExecuting)
  {
    throw ScanNotInitializedException();
  }

  if (this->attributeType == INTEGER){
    LeafNode<int>* currentNodeInt = (LeafNode<int>*) currentPageData;
    if(nextEntry == currentNodeInt->key_count || currentNodeInt->ridArray[nextEntry].page_number == 0)
    {
      // if invalid entry or at the end, try skip to right neighbor
      //if failed, then scan complete
      if(currentNodeInt->rightSibPageNo == 0)
      {
        throw IndexScanCompletedException();
      }
      currentPageNum = currentNodeInt->rightSibPageNo;
      this->bufMgr->readPage(file, currentPageNum, currentPageData);
      currentNodeInt = (LeafNode<int>*) currentPageData;
      this->bufMgr->unPinPage(file, currentPageNum, false);
      // nextEntry must be reset to 0
      nextEntry = 0;
    }
  
    // Check to make sure the key is in valid range
    int key = currentNodeInt->keyArray[nextEntry];

    if (lowOp == GTE && highOp == LTE && key >= this->lowValInt && key <= this->highValInt){
      outRid = currentNodeInt->ridArray[nextEntry];
      nextEntry++;
      return;
    } else if (lowOp == GTE && highOp == LT && key >= this->lowValInt && key < this->highValInt){
      outRid = currentNodeInt->ridArray[nextEntry];
      nextEntry++;
      return;
    } else if (lowOp == GT && highOp == LT && key > this->lowValInt && key < this->highValInt){
      outRid = currentNodeInt->ridArray[nextEntry];
      nextEntry++;
      return;
    } else if (lowOp == GT && highOp == LTE && key > this->lowValInt && key <= this->highValInt){
      outRid = currentNodeInt->ridArray[nextEntry];
      nextEntry++;
      return;
    } else {
      throw IndexScanCompletedException();
    }
      return;

  } else if (this->attributeType == DOUBLE){
    
      LeafNode<double>* currentNodeDouble = (LeafNode<double>*) currentPageData;
      if(nextEntry == currentNodeDouble->key_count || currentNodeDouble->ridArray[nextEntry].page_number == 0)
      {
        // if invalid entry or at the end, try skip to right neighbor
        //if failed, then scan complete
        if(currentNodeDouble->rightSibPageNo == 0)
        {
          throw IndexScanCompletedException();
        }
        currentPageNum = currentNodeDouble->rightSibPageNo;
        this->bufMgr->readPage(file, currentPageNum, currentPageData);
        currentNodeDouble = (LeafNode<double>*) currentPageData;
        this->bufMgr->unPinPage(file, currentPageNum, false);
        // nextEntry must be reset to 0
        nextEntry = 0;
      }
    
      // Check to make sure the key is in valid range
      double key = currentNodeDouble->keyArray[nextEntry];

      if (lowOp == GTE && highOp == LTE && key >= this->lowValDouble && key <= this->highValDouble){
        outRid = currentNodeDouble->ridArray[nextEntry];
        nextEntry++;
        return;
      } else if (lowOp == GTE && highOp == LT && key >= this->lowValDouble && key < this->highValDouble){
        outRid = currentNodeDouble->ridArray[nextEntry];
        nextEntry++;
        return;
      } else if (lowOp == GT && highOp == LT && key > this->lowValDouble && key < this->highValDouble){
        outRid = currentNodeDouble->ridArray[nextEntry];
        nextEntry++;
        return;
      } else if (lowOp == GT && highOp == LTE && key > this->lowValDouble && key <= this->highValDouble){
        outRid = currentNodeDouble->ridArray[nextEntry];
        nextEntry++;
        return;
      } else {
        throw IndexScanCompletedException();
      }

    } else {
        LeafNode<std::string>* currentNodeString = (LeafNode<std::string>*) currentPageData;
        int skipcount = 0;
        if(nextEntry == currentNodeString->key_count || currentNodeString->ridArray[nextEntry].page_number == 0)
        {
          // if invalid entry or at the end, try skip to right neighbor
          // if failed, then scan complete
          if(currentNodeString->rightSibPageNo == 0)
          {
            throw IndexScanCompletedException();
          }
          currentPageNum = currentNodeString->rightSibPageNo;
          this->bufMgr->readPage(file, currentPageNum, currentPageData);
          currentNodeString = (LeafNode<std::string>*) currentPageData;
          nextEntry = 0;
          this->bufMgr->unPinPage(file, currentPageNum, false);
        }
      
        // Check to make sure the key is in valid range
        std::string key;
        key = std::string(currentNodeString->keyArray[nextEntry]).substr(0, 10);
        // std::cout << "Compare key:  " << key <<"/  ";

        if (lowOp == GTE && highOp == LTE && (key >= this->lowValString) && (key <= this->highValString)){
          // std::cout << "M  ";
          outRid = currentNodeString->ridArray[nextEntry];
          nextEntry++;
          return;
        } else if (lowOp == GTE && highOp == LT && ((key >= this->lowValString)) && (key < this->highValString)){
          // std::cout << "M  " ;
          outRid = currentNodeString->ridArray[nextEntry];
          nextEntry++;
          return;
        } else if (lowOp == GT && highOp == LT && ((key > this->lowValString)) && (key < this->highValString)){
          // std::cout << "M  ";
          outRid = currentNodeString->ridArray[nextEntry];
          nextEntry++;
          return;
        } else if (lowOp == GT && highOp == LTE && ((key > this->lowValString)) && (key <= this->highValString)){
          // std::cout << "M  ";
          outRid = currentNodeString->ridArray[nextEntry];
          nextEntry++;
          return;
        } else {
          // std::cout << "SKIP COUNT:  " << skipcount <<"  ";
          throw IndexScanCompletedException();
        }
    }
    return;
}




// -----------------------------------------------------------------------------
// BTreeIndex::endScan
// -----------------------------------------------------------------------------

const void BTreeIndex::endScan() 
{
    // Ensure a scan is currently executing
    if (!scanExecuting) {
        throw ScanNotInitializedException();
    }

    // Reset scan state
    this->scanExecuting = false;
    
}

}
